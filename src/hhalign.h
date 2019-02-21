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
    HHalign(Parameters &par);

    virtual ~HHalign();

    virtual void run(FILE *query_fh, char *query_path);

    static void ProcessAllArguments(Parameters &par);

private:
    static void help(Parameters &par, char all = 0);

    static void ProcessArguments(Parameters &par);

    std::vector<std::string>& tfiles;
};

#endif /* HHALIGN_H_ */
