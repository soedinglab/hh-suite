/*
 * hhalign.h
 *
 *  Created on: Jun 24, 2014
 *      Author: meiermark
 */

#ifndef HHALIGN_H_
#define HHALIGN_H_

#include "hhblits.h"

class HHalign : public HHblits {
  public:
    HHalign(Parameters& par, std::vector<HHblitsDatabase*>& databases);
    virtual ~HHalign();
    void run(FILE* query_fh, char* query_path, std::vector<std::string>& template_paths);
    static void ProcessAllArguments(int argc, char** argv, Parameters& par);

  private:
    static void help(Parameters& par, char all=0);
    static void ProcessArguments(int argc, char** argv, Parameters& par);
};

#endif /* HHALIGN_H_ */
