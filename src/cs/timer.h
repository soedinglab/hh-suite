#ifndef CS_TIMER_H_
#define CS_TIMER_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctime>

namespace cs {

// Use time() call to keep track of elapsed time between creation and
// destruction.  If verbose is true, Timer will print a message showing
// elapsed time to the given output stream upon destruction.
// Adapted from timer class in bowtie source code.
class Timer {
  public:
    Timer(FILE* out = stdout, const char *msg = "", bool verbose = true) :
            t_(time(0)), out_(out), msg_(msg), verbose_(verbose) { }

    // Optionally print message
    ~Timer() {
        if(verbose_) Write(out_);
    }

    /// Return elapsed time since Timer object was created
    time_t elapsed() const {
        return time(0) - t_;
    }

    void Write(FILE* out) {
        time_t passed = elapsed();
        // Print the message supplied at construction time followed
        // by time elapsed formatted HH:MM:SS
        unsigned int hours   = (passed / 60) / 60;
        unsigned int minutes = (passed / 60) % 60;
        unsigned int seconds = (passed % 60);
        fprintf(out, "%s%02d:%02d:%02d\n", msg_, hours, minutes, seconds);
    }

  private:
    time_t      t_;
    FILE*       out_;
    const char* msg_;
    bool        verbose_;
};

static inline void LogTime(FILE* out, bool nl = true) {
    struct tm *current;
    time_t now;
    time(&now);
    current = localtime(&now);
    fprintf(out, "%02d:%02d:%02d", current->tm_hour, current->tm_min, current->tm_sec);
    if (nl) fputs("\n", out);
}

}  // namespace cs

#endif  // CS_TIMER_H_
