//
// Created by dtaliun on 1/31/20.
//

#include "aux.h"

//## some_file.geno[.gz] to some_file.site[.gz]
string build_sites_filename(const string& filename) {
    regex genogz_regex("\\.geno\\.gz$");
    regex geno_regex("\\.geno$");
    regex seqgz_regex("\\.seq\\.gz$");
    regex seq_regex("\\.seq$");
    if (regex_search(filename, genogz_regex)) {
        return regex_replace(filename, genogz_regex, ".site.gz");
    } else if (regex_search(filename, geno_regex)) {
        return regex_replace(filename, geno_regex, ".site");
    } else if (regex_search(filename, seqgz_regex)) {
        return regex_replace(filename, seqgz_regex, ".site.gz");
    } else if (regex_search(filename, seq_regex)) {
        return regex_replace(filename, seq_regex, ".site");
    } else {
        return filename + ".site";
    }
}

