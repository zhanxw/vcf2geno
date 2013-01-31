#include <cassert>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "tabix.h"

#include "Argument.h"
#include "Utils.h"
#include "VCFUtil.h"
#include "SimpleMatrix.h"
#include "IO.h"

/**
 * All indices are 0-based
 */
void appendGeno(SimpleMatrix* geno, const int peopleIndex, const int genoIndex, int g) {
  if (geno->nrow() <= peopleIndex) {
    geno->resize(geno->nrow() + 1, geno->ncol());
  }
  if (geno->ncol() <= genoIndex) {
    geno->resize(geno->nrow(), geno->ncol() + 1);
  }
  (*geno)[peopleIndex][genoIndex] = g;
}

int main(int argc, char** argv){
    time_t currentTime = time(0);
    fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

    ////////////////////////////////////////////////
    BEGIN_PARAMETER_LIST(pl)
        ADD_PARAMETER_GROUP(pl, "Input/Output")
        ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
        ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
        // ADD_PARAMETER_GROUP(pl, "People Filter")
        // ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "give IDs of people that will be included in study")
        // ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "from given file, set IDs of people that will be included in study")
        // ADD_STRING_PARAMETER(pl, peopleExcludeID, "--peopleExcludeID", "give IDs of people that will be included in study")
        // ADD_STRING_PARAMETER(pl, peopleExcludeFile, "--peopleExcludeFile", "from given file, set IDs of people that will be included in study")
        // ADD_PARAMETER_GROUP(pl, "Site Filter")
        // ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
        // ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
        END_PARAMETER_LIST(pl)
        ;    

    pl.Read(argc, argv);
    pl.Status();
    
    if (FLAG_REMAIN_ARG.size() > 0){
        fprintf(stderr, "Unparsed arguments: ");
        for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
            fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
        }
        fprintf(stderr, "\n");
        abort();
    }

    REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using: --inVcf");

    const char* fn = FLAG_inVcf.c_str(); 
    VCFInputFile vin(fn);

    // // set range filters here
    // // e.g.     
    // // vin.setRangeList("1:69500-69600");
    // vin.setRangeList(FLAG_rangeList.c_str());
    // vin.setRangeFile(FLAG_rangeFile.c_str());

    // // set people filters here
    // if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
    //     vin.excludeAllPeople();
    //     vin.includePeople(FLAG_peopleIncludeID.c_str());
    //     vin.includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
    // }
    // vin.excludePeople(FLAG_peopleExcludeID.c_str());
    // vin.excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());

    // store intemediate results
    OrderedMap < std::string, int> markerIndex;
    SimpleMatrix geno; // store genotypes by person
    
    // let's write it out.
    VCFOutputFile* vout = NULL;

    // if (FLAG_updateId != "") {
    //   int ret = vin.updateId(FLAG_updateId.c_str());
    //   fprintf(stdout, "%d samples have updated id.\n", ret);
    // }
    FILE* fSite = fopen( (FLAG_outPrefix + ".site").c_str(), "wt");
    fprintf(fSite, "CHROM\tPOS\tID\tREF\tALT\n");
    std::string markerName;
    int lineNo = 0;
    // int nonVariantSite = 0;
    while (vin.readRecord()){
        lineNo ++;
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        VCFIndividual* indv;
        // if (FLAG_variantOnly) {
        //   bool hasVariant = false;
        //   int geno;
        //   int GTidx = r.getFormatIndex("GT");
        //   for (int i = 0; i < people.size() ;i ++) {
        //     indv = people[i];
        //     geno = indv->justGet(0).getGenotype();
        //     if (geno != 0 && geno != MISSING_GENOTYPE)
        //       hasVariant = true;
        //   }
        //   if (!hasVariant) {
        //     nonVariantSite++;
        //     continue;
        //   }
        // }
        markerName = r.getID();
        if (markerName == ".") {
          markerName = r.getChrom();
          markerName += ":";
          markerName += r.getPosStr();
        }
        if (markerIndex.find(markerName)){
          fprintf(stderr, "Duplicated marker name [ %s ]\n", markerName.c_str());
          continue;
        }
        fprintf(fSite, "%s\t%s\t%s\t%s\t%s\n", r.getChrom(), r.getPosStr(), markerName.c_str(), r.getRef(), r.getAlt());

        const int index = markerIndex.size();
        markerIndex[markerName] = index;
            
        int GTidx = r.getFormatIndex("GT");
        bool missing;
        for (int i = 0; i < people.size(); ++i) {
          const VCFValue& v = people[i]->get(GTidx, &missing);
          if (!missing) {
            int g = v.getGenotype();
            appendGeno(&geno, i, index, g);
          } else {
            appendGeno(&geno, i, index, -9);
          }
        }
    };
    fclose(fSite);
    
    // output geno file
    FILE* fGeno = fopen ( (FLAG_outPrefix + ".geno").c_str(), "wt");
    std::vector<std::string> names;
    vin.getVCFHeader()->getPeopleName(&names);
    for (int i = 0; i < names.size(); ++i) {
      fprintf(fGeno, "%s\t%s", names[i].c_str(), names[i].c_str());
      for (int j = 0; j < geno.ncol(); ++j) {
        fprintf(fGeno, "\t%d", (int) geno[i][j]);
      }
      fprintf(fGeno, "\n");
    }
    fclose(fGeno);
    
    fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
    fprintf(stdout, "Total %d people and %d markers are outputted\n", geno.nrow(), geno.ncol());
    // if (FLAG_variantOnly) {
    //   fprintf(stdout, "Skipped %d non-variant VCF records\n", nonVariantSite);
    // }
    
    return 0; 
};
