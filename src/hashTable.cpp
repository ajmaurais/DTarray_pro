//
//  hashTable.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/7/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include <hashTable.hpp>

namespace hashTable{
	
	/**********************/
	/*     LinkedList     */
	/*********************/
	
	template<class _Tp>
	_Tp* LinkedList<_Tp>::consolidate(const _Tp& newItem)
	{
		for(typename std::list<_Tp>::iterator it = dat.begin();it != dat.end(); ++it)
		{
			if(*it == newItem)
			{
				(*it).consolidate(newItem);
				return &(*it);
			}
		}
		dat.push_back(newItem);
		return &(dat.back());
	}
	
	template<class _Tp>
	_Tp* LinkedList<_Tp>::getItem(std::string key)
	{
		for(typename std::list<_Tp>::iterator it = dat.begin();it != dat.end(); ++it)
		{
			if(*it == key)
				return &(*it);
		}
		return nullptr;
	}
	
	template<class _Tp>
	bool LinkedList<_Tp>::itemExists(std::string key) const
	{
		for(typename std::list<_Tp>::iterator it = dat.begin(); it != dat.end(); ++it)
			if(*it == key)
				return true;
		return false;
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::write(std::ofstream& outF, int fxnNum)
	{
		if(!outF)
			throw std::runtime_error("Error writing file. Bad std::ofstream!");
		
		if(dat.empty())
			return;
		
		for(typename std::list<_Tp>::iterator it = dat.begin();it != dat.end(); ++it)
			(*it).write(outF, fxnNum);
			
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::apply(int fxnNum)
	{
		if(dat.empty())
			return;
		
		for(typename std::list<_Tp>::iterator it = dat.begin(); it != dat.end(); ++it)
			(*it).apply(fxnNum);
	}
	
	/**********************/
	/*     HashTable     */
	/*********************/
	
	template<class _Tp>
	inline size_t HashTable<_Tp>::hash(std::string key) const
	{
		size_t h = FIRSTH;
		size_t sLen = key.length();
		for(size_t i = 0; i < sLen; i++)
			h = (h * A)^(key[i] * B);
		return h % size;
	}
	
	template<class _Tp>
	void HashTable<_Tp>::printHistogram() const
	{
		std::cout << std::endl << "Hash Table Contains ";
		std::cout << getLength() << " Items total" << std::endl;
		for (int i = 0; i < size; i++)
		{
			std::cout << i + 1 << ":\t";
			for(int j = 0; j < array[i].getLength(); j++)
				std::cout << " X";
			std::cout << std::endl;
		}
	}
	
	template<class _Tp>
	unsigned long HashTable<_Tp>::getLength() const
	{
		unsigned long c = 0;
		for (size_t i = 0; i < size; i++)
			c += array[i].getLength();
		
		return c;
	}
	
	template<class _Tp>
	void HashTable<_Tp>::write(std::ofstream& outF, int fxnNum)
	{
		if(!outF)
			throw std::runtime_error("Error writing file. Bad std::ofstream!");
		
		for(size_t i = 0; i < size; i++)
			array[i].write(outF, fxnNum);
	}
	
	template<class _Tp>
	void HashTable<_Tp>::apply(int fxnNum)
	{
		for(size_t i = 0; i < size; i++)
			array[i].apply(fxnNum);
	}
}//end of hashTable namespace

