#include "Analytics.h"
#include "Utility.h"
#include <limits>

struct ReportDiplotypeAssigment{
    double frequency;
    std::string reportID;
    std::string HT1;
    std::string HT2;
    std::string uniqueFreq;
};

class ReportDiplotypeAssigments {
    public:


    typedef std::unordered_map<std::string,ReportDiplotypeAssigment> ReportDiplotypeAssignmentList_t;
    typedef typename ReportDiplotypeAssignmentList_t::iterator iterator;
    typedef typename ReportDiplotypeAssignmentList_t::const_iterator c_iterator;

    explicit ReportDiplotypeAssigments() {};

    c_iterator begin() const {return ReportDiplotypeAssigmentList.cbegin();}
    c_iterator end() const {return ReportDiplotypeAssigmentList.cend();}

    bool update(const ReportDiplotypeAssigment & RDA) {
        iterator foundIt = ReportDiplotypeAssigmentList.find(RDA.reportID);
        if ( foundIt == ReportDiplotypeAssigmentList.end() ) {
            std::pair<std::string,ReportDiplotypeAssigment> insertee (RDA.reportID,RDA);
            insertee.second.uniqueFreq = "Yes";
            ReportDiplotypeAssigmentList.insert(insertee);
            return false;
        } else {
            if ( RDA.frequency - foundIt->second.frequency > ZERO) {
                foundIt->second = RDA;
                foundIt->second.uniqueFreq = "Yes";               
            } else if (abs(RDA.frequency - foundIt->second.frequency) < ZERO) {
                foundIt->second.uniqueFreq = "No";   
            }  
            return true;
        }
    }

    ReportDiplotypeAssigment find(const std::string & suche) {
        iterator foundIt = ReportDiplotypeAssigmentList.find(suche);
        return  foundIt->second;
    }

    private:

    ReportDiplotypeAssignmentList_t ReportDiplotypeAssigmentList; 
};

struct GenotypeDiplotypeAssignment{
    double frequency;
    std::string name;
    std::string HT1;
    std::string HT2;
};

class GenotypeDiplotypeAssignments {
    public:
    typedef std::unordered_map<size_t,GenotypeDiplotypeAssignment> GenotypeDiplotypeAssignmentList_t;
    typedef typename GenotypeDiplotypeAssignmentList_t::iterator iterator;
    typedef typename GenotypeDiplotypeAssignmentList_t::const_iterator c_iterator;

    explicit GenotypeDiplotypeAssignments() {};
    
    c_iterator c_listBegin() const {return genotypeDiplotypeAssignmentList.cbegin();}
    c_iterator c_listEnd() const {return genotypeDiplotypeAssignmentList.cend();}
        
    std::pair<iterator, bool> add(const GenotypeDiplotypeAssignment & gda){
        size_t hashValue = string_hash(gda.name);
        auto inserted = genotypeDiplotypeAssignmentList.emplace(hashValue, gda);
        return inserted;
    }
    
    std::pair<iterator, bool> addOrUpdate(const GenotypeDiplotypeAssignment & gda){
        size_t hashValue = string_hash(gda.name);
        auto inserted = genotypeDiplotypeAssignmentList.emplace(hashValue, gda);
        if (!inserted.second) {
            inserted.first->second.HT1 = gda.HT1;
            inserted.first->second.HT2 = gda.HT2;
            inserted.first->second.frequency = gda.frequency;
            inserted.second=true;
        }
        return inserted;
    }
    
    double findFrequencyByName (std::string searchName, double result = std::numeric_limits<double>::quiet_NaN()){
        for(auto iter = genotypeDiplotypeAssignmentList.cbegin();iter != genotypeDiplotypeAssignmentList.cend(); iter++){
            if (iter->second.name == searchName){
                result = iter->second.frequency;
                break;
            }
        }
        return result;
    }

private:
    GenotypeDiplotypeAssignmentList_t genotypeDiplotypeAssignmentList;
    std::hash<std::string> string_hash;
};

bool descending(std::pair<std::string,int> a, std::pair<std::string,int> b) {
    return a.second > b.second;
}

bool ascending(std::pair<std::string,int> a, std::pair<std::string,int> b) {
    return a.second < b.second;
}

void analyticsHaplotypeUsage( Phenotypes &pt, Haplotypes &ht , KeyPairs & kps, std::string  analyticsFileName){
    
    std::cout << "#########Running Analytics" << std::endl;

    const std::string tab = "\t";
    
    GenotypeDiplotypeAssignments probableDiplotypes;
    GenotypeDiplotypeAssignment  genotypeDiplotype;
    ReportDiplotypeAssigments reports; 
    ReportDiplotypeAssigment report;
    
    std::unordered_map<double,std::string> freqPerReportList;
    std::unordered_map<std::string, double> freqSumPerReport;
    std::unordered_map<std::string, int> haplotypeUseCount;
    std::hash<std::string> string_hash;
    
    for (auto itPheno = pt.c_listBegin();itPheno != pt.c_listEnd();itPheno++) {
        
        string_vector_t reportIDList = itPheno->second.getReportID_all();

        double maxFrequency = 0.0;
        std::string maxHT1 = "N/Found";
        std::string  maxHT2 = "N/Found";
        
        for(auto itDiplo = itPheno->second.c_diplotypeListBegin();itDiplo != itPheno->second.c_diplotypeListEnd(); itDiplo ++) {
            
            for (auto reportID : reportIDList) {freqSumPerReport[reportID] += itDiplo->frequency;}
            
            if ( itDiplo->frequency > maxFrequency ) {
                 maxFrequency = itDiplo->frequency;
                 maxHT1 = kps.findNameByID(itDiplo->id1);
                 maxHT2 = kps.findNameByID(itDiplo->id2);
            }
        }

        for (auto reportID : reportIDList) {
            if (probableDiplotypes.findFrequencyByName(reportID,0) < maxFrequency ) {
                genotypeDiplotype.frequency  = maxFrequency;
                genotypeDiplotype.name = reportID;
                genotypeDiplotype.HT1 = maxHT1;
                genotypeDiplotype.HT2 = maxHT2;
                probableDiplotypes.addOrUpdate(genotypeDiplotype);

                report.frequency = maxFrequency;
                report.reportID = reportID;
                report.HT1 = maxHT1;
                report.HT2 = maxHT2;
                reports.update(report);
            }
        }
    }
    
    for (auto itA = probableDiplotypes.c_listBegin()   ;itA!=probableDiplotypes.c_listEnd();itA++){
        haplotypeUseCount[itA->second.HT1]++;
        haplotypeUseCount[itA->second.HT2]++;
    }
    
    std::vector<std::pair<std::string,int>> sortedVector(haplotypeUseCount.begin(),haplotypeUseCount.end());
    std::sort(sortedVector.begin(), sortedVector.end(), descending);
    
    int intSum=0;
    for (auto hUC : sortedVector) intSum+=hUC.second;
    
    std::ofstream outFile;
    openFileToWrite(analyticsFileName+"_1.txt",outFile);
    
    for (auto hUC : sortedVector){
        outFile << hUC.first << tab
                  << hUC.second << tab
                  << hUC.second*1./intSum << tab
                  << ht.getFrequency(string_hash(hUC.first))
                  << std::endl ;
    }
    outFile.close();

    std::vector<std::string> sortedReportIDs;
    for (auto iterat : reports) {
        sortedReportIDs.push_back(iterat.first);
    }
       
    std::sort(sortedReportIDs.begin(),sortedReportIDs.end(),[](const std::string& a,const std::string& b) {
        if (a.size() == b.size()) {
            return a < b;            
        } else if (a.size() > b.size()) {
            return false;
        } else {
            return true;
        }
    });

    openFileToWrite(analyticsFileName+"_2.txt",outFile);
    for (std::string repID : sortedReportIDs) {
        auto rep = reports.find(repID);
            outFile << rep.reportID << tab
                    << rep.HT1 << tab
                    << rep.HT2 << tab
                    << rep.uniqueFreq 
                    << std::endl;
    }
    outFile.close();
}
