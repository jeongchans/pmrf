#ifndef _SEQIO_H_
#define _SEQIO_H_

#include <istream>
#include <string>

class SeqRecord {
  public:
    SeqRecord() {};
    SeqRecord(const std::string& id, const std::string& desc, const std::string& seq) : id(id), desc(desc), seq(seq) {};

    std::string id;
    std::string desc;
    std::string seq;
};

class SeqParser {
  public:
    virtual SeqRecord next() = 0;
    virtual bool has_next() = 0;
};

class FastaParser : public SeqParser {
  public:
    FastaParser(std::istream& is) : is(is) { init(); }

    virtual SeqRecord next();
    virtual bool has_next() { return _has_next; }

  private:
    std::istream& is;   // input stream
    bool _has_next;
    std::string buffer;

    void init();
    inline void set_buffer(const std::string& s) { buffer = s; _has_next = true; }
};

#endif
