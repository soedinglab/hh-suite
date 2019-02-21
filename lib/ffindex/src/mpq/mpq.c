#include "mpq.h"

int MPQ_size;
int MPQ_rank;
size_t MPQ_num_jobs;
int MPQ_is_init = 0;

enum {
    TAG_JOB,
    TAG_FREE,
    TAG_FINISHED
};

enum {
    MSG_RELEASE,
    MSG_JOB,
    MSG_FREE,
    MSG_FINISHED
};

int MPQ_Init(int argc, const char **argv, const size_t num_jobs) {
    if (MPQ_is_init == 1) {
        return MPQ_ERROR_REINIT;
    }

    MPI_Init(&argc, (char***)&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPQ_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPQ_size);

    MPQ_num_jobs = num_jobs;

    int workers = MPQ_size - 1;
    if (workers < 1) {
        return MPQ_ERROR_NO_WORKERS;
    }

    MPQ_is_init = 1;

    return MPQ_SUCCESS;
}

void MPQ_Worker(MPQ_Payload_t payload, void *env) {
    if (MPQ_is_init == 0) {
        return;
    }

    int message_free = MSG_FREE;
    int message_finished = MSG_FINISHED;
    int message_job[3];

    while (1) {
        MPI_Send(&message_free, 1, MPI_INT, MPQ_MASTER, TAG_FREE, MPI_COMM_WORLD);

        MPI_Recv(message_job, 3, MPI_INT, MPQ_MASTER, TAG_JOB, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (message_job[0] == MSG_RELEASE) {
            MPI_Send(&message_finished, 1, MPI_INT, MPQ_MASTER, TAG_FINISHED, MPI_COMM_WORLD);
            break;
        }

        payload(env, message_job[1], message_job[2]);
    }
}

void MPQ_Worker_Release(const int worker_rank) {
    int message[3];
    message[0] = MSG_RELEASE;
    message[1] = 0;
    message[2] = 0;
    MPI_Send(message, 3, MPI_INT, worker_rank, TAG_JOB, MPI_COMM_WORLD);
}

void MPQ_Worker_Job(const int worker_rank, int start, int end) {
    int message[3];
    message[0] = MSG_JOB;
    message[1] = start;
    message[2] = end;
    MPI_Send(message, 3, MPI_INT, worker_rank, TAG_JOB, MPI_COMM_WORLD);
}

void MPQ_Master(const size_t split_size) {
    if (MPQ_is_init == 0) {
        return;
    }

    size_t split = split_size;
    if (MPQ_num_jobs < split) {
        split = 1;
    }

    MPI_Status status;
    int message[3];

    for (size_t i = 0; i < MPQ_num_jobs; i += split) {
        MPI_Recv(message, 1, MPI_INT, MPI_ANY_SOURCE, TAG_FREE, MPI_COMM_WORLD, &status);
        int worker_rank = status.MPI_SOURCE;

        size_t end = (MPQ_num_jobs < (i + split)) ? MPQ_num_jobs : (i + split);
        MPQ_Worker_Job(worker_rank, i, end);
    }

    for (int i = 1; i < MPQ_size; i++) {
        MPQ_Worker_Release(i);
        MPI_Recv(message, 1, MPI_INT, i, TAG_FINISHED, MPI_COMM_WORLD, &status);
    }
}
