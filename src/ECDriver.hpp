/***
 **
 *  File: ECDriver.h
 *  Created: Dec 12, 2009
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *
 *  Copyright 2009 Xiao Yang, Karin Dorman, Srinivas Aluru
 *  Copyright 2011-12 Ankit Shah, Sriram P C, Srinivas Aluru
 *
 *  This file is part of Reptile (version 1.1).
 *
 *  Reptile is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Reptile is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _ECDRIVER_H
#define	_ECDRIVER_H


#include "util.h"
#include "ECData.hpp"

#include <fstream>
#include <string>

class ECDriver {
public:
    ECDriver(ECData& p, Para& mpara):
        ecdata_(p),inPara_(mpara),readID_(0){};
    virtual ~ECDriver(){};

    void ec();
private:
    ECData& ecdata_;
    Para& inPara_;

    int readID_;
    std::ofstream oHandle_;

    void correctBatch(const ReadStore& rbatch);

    void correctBatchST(const ReadStore& rbatch);

    void correctBatchMT(const ReadStore& rbatch);

    void correctReadsFromFile();

    void ecStatic();

    void ecDynamic();
};

#endif	/* _ECDRIVER_H */
