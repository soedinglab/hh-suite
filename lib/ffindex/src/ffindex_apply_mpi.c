/* 
 * FFindex
 * written by Andy Hauser <hauser@genzentrum.lmu.de>,
 * and Milot Mirdita <milot@mirdita.de>
 * Please add your name here if you distribute modified versions.
 *
 * FFindex is provided under the Create Commons license "Attribution-ShareAlike
 * 3.0", which basically captures the spirit of the Gnu Public License (GPL).
 *
 * See:
 * http://creativecommons.org/licenses/by-sa/3.0/
 *
 * ffindex_apply_mpi
 * apply a program to each FFindex entry
 */

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <limits.h> // PIPE_BUF
#include <stdlib.h> // EXIT_*, system, malloc, free
#include <unistd.h> // pipe, fork, close, dup2, execvp, write, read, opt*

#include <sys/time.h>
#include <sys/wait.h>
#include <sys/mman.h> // munmap
#include <fcntl.h>    // fcntl, F_*, O_*
#include <signal.h>   // sigaction, sigemptyset

#include <getopt.h>   // getopt_long

#include <spawn.h>     // spawn_*

#include "ffindex.h"
#include "ffutil.h"

#ifdef HAVE_MPI
#include "mpq/mpq.h"
#endif

char read_buffer[400 * 1024 * 1024];
extern char **environ;

int
ffindex_apply_by_entry(char *data, ffindex_entry_t *entry, char *program_name, char **program_argv,
                       FILE *data_file_out, FILE *index_file_out, FILE *log_file_out, size_t *offset, int quiet) {
    int ret = 0;
    int capture_stdout = (data_file_out != NULL && index_file_out != NULL);

    int pipefd_stdin[2];
    int pipefd_stdout[2];

    ret = pipe(pipefd_stdin);
    if (ret != 0) {
        fprintf(stderr, "ERROR in pipe stdin!\n");
        perror(entry->name);
        return errno;
    }

    if (capture_stdout) {
        ret = pipe(pipefd_stdout);
        if (ret != 0) {
            fprintf(stderr, "ERROR in pipe stdout!\n");
            perror(entry->name);
            return errno;
        }

        // Flush so child doesn't copy and also flushes, leading to duplicate output
        fflush(data_file_out);
        fflush(index_file_out);
    }

    posix_spawn_file_actions_t factions;
    posix_spawn_file_actions_init(&factions);
    posix_spawn_file_actions_addclose(&factions, pipefd_stdin[1]);
    posix_spawn_file_actions_adddup2(&factions, pipefd_stdin[0], fileno(stdin));
    if (capture_stdout) {
        posix_spawn_file_actions_addclose(&factions, pipefd_stdout[0]);
        posix_spawn_file_actions_adddup2(&factions, pipefd_stdout[1], fileno(stdout));
    }

    short spawnFlags = 0;
    posix_spawnattr_t attr;
    posix_spawnattr_init(&attr);
#ifdef POSIX_SPAWN_USEVFORK
    spawnFlags |= POSIX_SPAWN_USEVFORK;
#endif
    posix_spawnattr_setflags(&attr, spawnFlags);

    pid_t child_pid;

    struct timeval tv;
    gettimeofday(&tv, NULL);
    int64_t start = (tv.tv_sec) * 1000LL + (tv.tv_usec) / 1000;
    int err = posix_spawnp(&child_pid, program_name, &factions, &attr, program_argv, environ);
    if (err) {
        fprintf(stderr, "ERROR in fork()\n");
        perror(entry->name);
        return errno;
    }
    else {
        // parent writes to and possible reads from child
        int flags = 0;

        // Read end is for child only
        close(pipefd_stdin[0]);

        if (capture_stdout) {
            close(pipefd_stdout[1]);
        }

        char *file_data = ffindex_get_data_by_entry(data, entry);

        if (capture_stdout) {
            flags = fcntl(pipefd_stdout[0], F_GETFL, 0);
            fcntl(pipefd_stdout[0], F_SETFL, flags | O_NONBLOCK);
        }

        char *b = read_buffer;

        // Write file data to child's stdin.
        ssize_t written = 0;
        // Don't write ffindex trailing '\0'
        size_t to_write = entry->length - 1;
        while (written < to_write) {
            size_t rest = to_write - written;
            size_t batch_size = PIPE_BUF;
            if (rest < PIPE_BUF) {
                batch_size = rest;
            }

            ssize_t w = write(pipefd_stdin[1], file_data + written, batch_size);
            if (w < 0 && errno != EPIPE) {
                fprintf(stderr, "ERROR in child!\n");
                perror(entry->name);
                break;
            }
            else {
                written += w;
            }

            if (capture_stdout) {
                // To avoid blocking try to read some data
                ssize_t r = read(pipefd_stdout[0], b, PIPE_BUF);
                if (r > 0) {
                    b += r;
                }
            }
        }
        close(pipefd_stdin[1]); // child gets EOF

        if (capture_stdout) {
            // Read rest
            fcntl(pipefd_stdout[0], F_SETFL, flags); // Remove O_NONBLOCK
            ssize_t r;
            while ((r = read(pipefd_stdout[0], b, PIPE_BUF)) > 0) {
                b += r;
            }
            close(pipefd_stdout[0]);

            ffindex_insert_memory(data_file_out, index_file_out, offset, read_buffer, b - read_buffer, entry->name);
        }

        int status;
        waitpid(child_pid, &status, 0);
        if (!quiet) {
            struct timeval tv;
            gettimeofday(&tv, NULL);
            int64_t end = (tv.tv_sec) * 1000LL + (tv.tv_usec) / 1000;
            fprintf(log_file_out, "%s\t%ld\t%ld\t%ld\t%d\n", entry->name, entry->offset, entry->length, end - start, WEXITSTATUS(status));
        }
    }

    return EXIT_SUCCESS;
}

#ifdef HAVE_MPI
typedef struct ffindex_apply_mpi_data_s ffindex_apply_mpi_data_t;
struct ffindex_apply_mpi_data_s {
    void *data;
    ffindex_index_t *index;
    FILE *data_file_out;
    FILE *index_file_out;
    FILE *log_file_out;
    char *program_name;
    char **program_argv;
    int quiet;
    size_t offset;
};

void ffindex_apply_worker_payload(void *pEnv, const size_t start, const size_t end) {
    ffindex_apply_mpi_data_t *env = (ffindex_apply_mpi_data_t *) pEnv;
    for (size_t i = start; i < end; i++) {
        ffindex_entry_t *entry = ffindex_get_entry_by_index(env->index, i);
        if (entry == NULL) {
            break;
        }

        int error = ffindex_apply_by_entry(env->data, entry,
                                           env->program_name,
                                           env->program_argv,
                                           env->data_file_out,
                                           env->index_file_out,
                                           env->log_file_out,
                                           &(env->offset),
                                           env->quiet);
        if (error != 0) {
            perror(entry->name);
            break;
        }
    }

    return;
}

#endif

void ignore_signal(int signal) {
    struct sigaction handler;
    handler.sa_handler = SIG_IGN;
    sigemptyset(&handler.sa_mask);
    handler.sa_flags = 0;
    sigaction(signal, &handler, NULL);
}

void usage() {
    fprintf(stderr,
            "USAGE: ffindex_apply_mpi [-q] [-k] "
#ifdef HAVE_MPI
                    "[-p PARTS] [-l LOG_FILENAME_PREFIX] "
#endif
                    "[-d DATA_FILENAME_OUT -i INDEX_FILENAME_OUT] DATA_FILENAME INDEX_FILENAME -- PROGRAM [PROGRAM_ARGS]*\n"
                    "\nDesigned and implemented by Andy Hauser <hauser@genzentrum.lmu.de> and Milot Mirdita <milot@mirdita.de>.\n\n"
#ifdef HAVE_MPI
                    "\t[-p PARTS]\t\tSets how many entries one worker processes per job.\n"
                    "\t[-l LOG_FILE_PREFIX]\tPrefix for filename for the per worker process logfiles.\n"
#endif
                    "\t[-q]\t\t\tSilence the logging of every processed entry.\n"
                    "\t[-k]\t\t\tKeep unmerged ffindex splits.\n"
                    "\t[-d DATA_FILENAME_OUT]\tFFindex data file where the results will be saved to.\n"
                    "\t[-i INDEX_FILENAME_OUT]\tFFindex index file where the results will be saved to.\n"
                    "\tDATA_FILENAME\t\tInput ffindex data file.\n"
                    "\tINDEX_FILENAME\t\tInput ffindex index file.\n"
                    "\tPROGRAM [PROGRAM_ARGS]\tProgram to be executed for every ffindex entry.\n"
    );
}

int main(int argn, char **argv) {
    ignore_signal(SIGPIPE);

    int exit_status = EXIT_SUCCESS;

    int quiet = 0;
    int keepTmp = 0;
    char *data_filename_out = NULL;
    char *index_filename_out = NULL;

#ifdef HAVE_MPI
    size_t parts = 1;
    char *log_filename = NULL;
#endif

    static struct option long_options[] =
            {
#ifdef HAVE_MPI
                    {"parts", required_argument, NULL, 'p'},
                    {"logfile", required_argument, NULL, 'l'},
#endif
                    {"data", required_argument, NULL, 'd'},
                    {"index", required_argument, NULL, 'i'},
                    {"quiet", no_argument, NULL, 'q'},
                    {"keep-tmp", no_argument, NULL, 'k'},
                    {NULL, 0, NULL, 0}
            };

    int opt;
    while (1) {
        int option_index = 0;
#ifdef HAVE_MPI
        const char *short_options = "kql:p:d:i:";
#else
        const char* short_options = "kqd:i:";
#endif
        opt = getopt_long(argn, argv, short_options, long_options, &option_index);

        if (opt == -1)
            break;

        switch (opt) {
#ifdef HAVE_MPI
            case 'p':
                parts = strtoull(optarg, NULL, 10);
                break;
            case 'l':
                log_filename = optarg;
                break;
#endif
            case 'd':
                data_filename_out = optarg;
                break;
            case 'i':
                index_filename_out = optarg;
                break;
            case 'q':
                quiet = 1;
                break;
            case 'k':
                keepTmp = 1;
                break;
            default:
                break;
        }
    }

    const int remaining_arguments = argn - optind;
    if (remaining_arguments < 3) {
        if (remaining_arguments < 2) {
            fprintf(stderr, "Please specify input data and index file.\n\n");
        }

        fprintf(stderr, "Please specify a program to execute.\n\n");

        usage();

        return EXIT_FAILURE;
    }

    if ((!data_filename_out && index_filename_out) || (data_filename_out && !index_filename_out)) {
        fprintf(stderr, "Please specify both output data and index file.\n\n");

        usage();

        return EXIT_FAILURE;
    }

    char *data_filename = argv[optind++];
    FILE *data_file = fopen(data_filename, "r");
    if (data_file == NULL) {
        fferror_print(__FILE__, __LINE__, argv[0], data_filename);
        return EXIT_FAILURE;
    }

    char *index_filename = argv[optind++];
    FILE *index_file = fopen(index_filename, "r");
    if (index_file == NULL) {
        fferror_print(__FILE__, __LINE__, argv[0], index_filename);
        exit_status = EXIT_FAILURE;
        goto cleanup_1;
    }

    char *program_name = argv[optind];
    char **program_argv = argv + optind;

    size_t data_size;
    char *data = ffindex_mmap_data(data_file, &data_size);
    if (data == MAP_FAILED) {
        fferror_print(__FILE__, __LINE__, "ffindex_mmap_data", data_filename);
        exit_status = EXIT_FAILURE;
        goto cleanup_2;
    }

    size_t entries = ffcount_lines(index_filename);
    ffindex_index_t *index = ffindex_index_parse(index_file, entries);
    if (index == NULL) {
        fferror_print(__FILE__, __LINE__, "ffindex_index_parse", index_filename);
        exit_status = EXIT_FAILURE;
        goto cleanup_3;
    }

#ifdef HAVE_MPI
    if (log_filename != NULL && quiet) {
        fprintf(stderr, "Please specify either quiet (-q) or a logfile (-l).\n\n");
        usage();
        exit_status = EXIT_FAILURE;
        goto cleanup_3;
    }

    int mpq_status = MPQ_Init(argn, argv, index->n_entries);
    if (mpq_status == MPQ_SUCCESS) {
        if (MPQ_rank != MPQ_MASTER) {
            ffindex_apply_mpi_data_t *env = malloc(sizeof(ffindex_apply_mpi_data_t));
            env->data = data;
            env->index = index;
            env->program_name = program_name;
            env->program_argv = program_argv;
            env->quiet = quiet;
            env->offset = 0;

            env->data_file_out = NULL;
            if (data_filename_out != NULL) {
                char data_filename_out_rank[FILENAME_MAX];
                snprintf(data_filename_out_rank, FILENAME_MAX, "%s.%d",
                         data_filename_out, MPQ_rank);

                env->data_file_out = fopen(data_filename_out_rank, "w+");
                if (env->data_file_out == NULL) {
                    fferror_print(__FILE__, __LINE__, "fopen", data_filename_out_rank);
                    exit_status = EXIT_FAILURE;
                    goto cleanup_4;
                }
            }

            env->index_file_out = NULL;
            if (index_filename_out != NULL) {
                char index_filename_out_rank[FILENAME_MAX];
                snprintf(index_filename_out_rank, FILENAME_MAX, "%s.%d",
                         index_filename_out, MPQ_rank);

                env->index_file_out = fopen(index_filename_out_rank, "w+");
                if (env->index_file_out == NULL) {
                    fferror_print(__FILE__, __LINE__, "fopen", index_filename_out_rank);
                    exit_status = EXIT_FAILURE;
                    goto cleanup_5;
                }
            }

            env->log_file_out = stdout;
            if (log_filename != NULL) {
                char log_filename_out_rank[FILENAME_MAX];
                snprintf(log_filename_out_rank, FILENAME_MAX, "%s.%d",
                         log_filename, MPQ_rank);

                env->log_file_out = fopen(log_filename_out_rank, "w+");
                if (env->log_file_out == NULL) {
                    fferror_print(__FILE__, __LINE__, "fopen", log_filename_out_rank);
                    exit_status = EXIT_FAILURE;
                    goto cleanup_6;
                }
            }

            MPQ_Worker(ffindex_apply_worker_payload, env);

            // make sure that all written files are properly flushed and synced
            // so that ffmerge_splits wont work on stale data
            if (env->log_file_out) {
                int fd = fileno(env->log_file_out);
                fflush(env->log_file_out);
                fsync(fd);
                fclose(env->log_file_out);
            }
            cleanup_6:
            if (env->index_file_out) {
                int fd = fileno(env->index_file_out);
                fflush(env->index_file_out);
                fsync(fd);
                fclose(env->index_file_out);
            }
            cleanup_5:
            if (env->data_file_out) {
                int fd = fileno(env->data_file_out);
                fflush(env->data_file_out);
                fsync(fd);
                fclose(env->data_file_out);
            }
            cleanup_4:
            free(env);

            MPI_Barrier(MPI_COMM_WORLD);
        } else {
            MPQ_Master(parts);
            MPI_Barrier(MPI_COMM_WORLD);

            int removeTmp = keepTmp == 0;
            ffmerge_splits(data_filename_out, index_filename_out, 1, MPQ_size - 1, removeTmp);
        }
    } else {
        if (mpq_status == MPQ_ERROR_NO_WORKERS) {
            fprintf(stderr, "MPQ_Init: Needs at least one worker process.\n");
            exit_status = EXIT_FAILURE;
        }
    }
#else
    FILE* data_file_out = NULL;
    if (data_filename_out != NULL)
    {
        char data_filename_out_rank[FILENAME_MAX];
        snprintf(data_filename_out_rank, FILENAME_MAX, "%s", data_filename_out);

        data_file_out = fopen(data_filename_out_rank, "w+");
        if (data_file_out == NULL)
        {
            fferror_print(__FILE__, __LINE__, "ffindex_apply_worker_payload", data_filename_out);
            return EXIT_FAILURE;
        }
    }

    FILE* index_file_out = NULL;
    if (index_filename_out != NULL)
    {
        char index_filename_out_rank[FILENAME_MAX];
        snprintf(index_filename_out_rank, FILENAME_MAX, "%s", index_filename_out);

        index_file_out = fopen(index_filename_out_rank, "w+");
        if (index_file_out == NULL)
        {
            fferror_print(__FILE__, __LINE__, "ffindex_apply_worker_payload", index_filename_out);
            return EXIT_FAILURE;
        }
    }

    size_t offset = 0;
    for (size_t i = 0; i < index->n_entries; i++)
    {
        ffindex_entry_t *entry = ffindex_get_entry_by_index(index, i);
        if (entry == NULL)
        {
            exit_status = errno;
            break;
        }

        int error = ffindex_apply_by_entry(data, entry,
                                           program_name, program_argv,
                                           data_file_out, index_file_out, stdout,
                                           &offset, quiet);
        if (error != 0)
        {
            perror(entry->name);
            exit_status = errno;
            break;
        }
    }

    if (index_file_out) {
        fclose(index_file_out);
    }
    if (data_file_out) {
        fclose(data_file_out);
    }
#endif

    munmap(index->index_data, index->index_data_size);
    free(index);

    cleanup_3:
    munmap(data, data_size);

    cleanup_2:
    if (index_file) {
        fclose(index_file);
    }

    cleanup_1:
    if (data_file) {
        fclose(data_file);
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return exit_status;
}
