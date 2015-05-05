#ifndef ALIGNMENT_HPP_
#define ALIGNMENT_HPP_

#include <string>

#include "sequence.hpp"

namespace aln
{
  class Alignment
  {
  private:
    aln::Sequence sequence1_;
    aln::Sequence sequence2_;
  public:
    Alignment(aln::Sequence& sequence2, aln::Sequence& sequence2):
      sequence1_(sequence1),
      sequence2_(sequence2){}
    
    aln::Sequence sequence1(){ return this->sequence1_; };
    aln::Sequence sequence2(){ return this->sequence2_; };
  }
}

#endif /* ALIGNMENT_HPP_ */