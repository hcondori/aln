#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include <string>

namespace aln
{
  class Sequence
  {
  private:
    std::string id_;
    std::string data_;
  public:
    Sequence(std::string id, std::string data):
      id_(id),
      data_(data){}
    Sequence(char* id, char* data):
      id_(id),
      data_(data){}
    std::string id(){ return this->id_; };
    std::string data(){ return this->data_; };
  };
}

#endif /* SEQUENCE_HPP_ */