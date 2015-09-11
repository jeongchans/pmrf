#ifndef _PROTBINFO_HMMIO_H_
#define _PROTBINFO_HMMIO_H_

#include <istream>
#include <gtest/gtest_prod.h>

#include "profilehmm.h"

//enum ProfileHMMFormat { HHM };

class HMMExporter {
  public:
    virtual void export_model(const ProfileHMM& model, std::ostream& os) = 0;
};

class HHsearchHMMExporter : public HMMExporter {
  public:
    HHsearchHMMExporter() {};

    virtual void export_model(const ProfileHMM& model, std::ostream& os);

  protected:
    inline void export_elem(const double& x, std::ostream& os);
    inline void export_elem(const char& x, std::ostream& os);

    FRIEND_TEST(HHsearchHMMExporterTest, test_export_emit_symbol);
    FRIEND_TEST(HHsearchHMMExporterTest, test_export_emit);
    FRIEND_TEST(HHsearchHMMExporterTest, test_export_transit);
    FRIEND_TEST(HHsearchHMMExporterTest, test_export_eff_num);

    void export_header(const ProfileHMM& model, std::ostream& os);
    void export_body(const ProfileHMM& model, std::ostream& os);
    void export_seq(const std::string& seq, const size_t& width, std::ostream& os);
    void export_emit_symbol(const std::string& sym, std::ostream& os);
    void export_emit(const Float1dArray& v, std::ostream& os);
    void export_transit(const StateType& type, const Float1dArray& v, std::ostream& os);
    void export_eff_num(const double& v, std::ostream& os);
};

class HMMImporter {
  public:
    virtual ProfileHMM import_model(std::istream& is) = 0;
};

class HHsearchHMMImporter : public HMMImporter {
  public:
    HHsearchHMMImporter(const Alphabet& abc) : abc(abc) {};

    virtual ProfileHMM import_model(std::istream& is);

  protected:
    const Alphabet& abc;

    inline double import_elem(std::istream& is);

    FRIEND_TEST(HHsearchHMMImporterTest, test_import_emit_symbol);
    FRIEND_TEST(HHsearchHMMImporterTest, test_import_emit);
    FRIEND_TEST(HHsearchHMMImporterTest, test_import_transit);
    FRIEND_TEST(HHsearchHMMImporterTest, test_import_eff_num);

    ProfileHMM import_header(std::istream& is);
    void import_body(ProfileHMM& model, std::istream& is);
    std::string import_seq(const size_t& length, std::istream& is);
    std::string import_emit_symbol(std::istream& is);
    Float1dArray import_emit(std::istream& is, const int& size);
    Float1dArray import_transit(std::istream& is, const StateType& type);
    double import_eff_num(std::istream& is);
};

//class NoEntryException : public std::exception {
//  public:
//    virtual const char* what() const throw() {
//        return "No HCI model found";
//    }
//};

//class HCIModelImporter {
//  public:
//    HCIModelImporter(const Alphabet& abc) : abc(abc) {};
//
//    HCIModel import_model(std::istream& is);
//
//  private:
//    const Alphabet& abc;
//    inline double import_elem(std::istream& is);
//
//    FRIEND_TEST(HCIModelImporterTest, test_import_emit_symbol);
//    FRIEND_TEST(HCIModelImporterTest, test_import_emit);
//    FRIEND_TEST(HCIModelImporterTest, test_import_transit);
//    FRIEND_TEST(HCIModelImporterTest, test_import_eff_num);
//    FRIEND_TEST(HCIModelImporterTest, test_extract_hcisvm_content);
//
//    std::string import_emit_symbol(std::istream& is);
//    Float1dArray import_emit(std::istream& is, const int& size);
//    Float1dArray import_transit(std::istream& is, const StateType& type);
//    double import_eff_num(std::istream& is);
//    Float2dArray import_cm_avg(std::istream& is, const int& size);
//    string extract_hcisvm_content(std::istream& is);
//};

#endif
