/*
 Copyright (C) 2015 Héctor Condori Alagón.

 This file is part of ALN, the massive Smith-Waterman pairwise sequence aligner.

 ALN is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ALN is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ALN.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ALN_BUFFER_HPP_
#define ALN_BUFFER_HPP_

namespace aln
{
  
  template<typename T>
  class Buffer
  {
  private:
    T* data_;
    int size_;
  public:
    Buffer(size_t size);
    Buffer(size_t size, size_t alignment);
    Buffer(const Buffer& buffer);
    ~Buffer();
    T* data();
    void resize(int size);
    T& operator [](const int index);
    Buffer& operator=(const Buffer& buffer);
  };

  template<typename T>
  Buffer<T>::Buffer(size_t size)
  {
    this->size_ = size;
    this->data_ = (T*)malloc(size * sizeof(T));
  }
  
  template<typename T>
  Buffer<T>::Buffer(size_t size, size_t alignment)
  {
    this->size_ = size;
    this->data_ = (T*)aligned_alloc(alignment, size * sizeof(T));
  }
  
  template<typename T>
  Buffer<T>::Buffer(const Buffer<T>& buffer)
  {
    this->size_ = buffer.size_;
    this->data_ = buffer.data_;
    //buffer.data_ = nullptr;
  }
  
  template<typename T>
  Buffer<T>& Buffer<T>::operator=(const Buffer<T>& buffer)
  {
    this-> data_ = buffer.data_;
    this-> size_ = buffer.size_;
    buffer.data_ = nullptr;
    return *this;
  }

  template<typename T>
  Buffer<T>::~Buffer()
  {
    free(this->data_);
  }

  template<typename T>
  T* Buffer<T>::data()
  {
    return this->data_;
  }

  template<typename T>
  T& Buffer<T>::operator [](const int index)
  {
    return this->data_[index];
  }
  
}

#endif /* ALN_BUFFER_HPP_ */