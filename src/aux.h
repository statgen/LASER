//
// Created by dtaliun on 1/31/20.
//

#ifndef LASER_AUX_H
#define LASER_AUX_H

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <regex>
#include <zlib.h>

using namespace std;

//## some_file.geno[.gz] to some_file.site[.gz]
string build_sites_filename(const string& filename);

#endif //LASER_AUX_H
