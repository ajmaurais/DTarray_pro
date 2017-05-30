//
//  hashTable.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/7/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "hashTable.hpp"

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
	_Tp* LinkedList<_Tp>::getItem(string key)
	{
		for(typename list<_Tp>::iterator it = dat.begin();it != dat.end(); ++it)
		{
			if(*it == key)
				return &(*it);
		}
		return nullptr;
	}
	
	template<class _Tp>
	bool LinkedList<_Tp>::itemExists(string key) const
	{
		for(typename list<_Tp>::iterator it = dat.begin(); it != dat.end(); ++it)
			if(*it == key)
				return true;
		return false;
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::write(ofstream& outF, int fxnNum)
	{
		if(!outF)
			throw runtime_error("Error writing file. Bad ofstream!");
		
		if(dat.empty())
			return;
		
		for(typename list<_Tp>::iterator it = dat.begin();it != dat.end(); ++it)
			(*it).write(outF, fxnNum);
			
	}
	
	/**********************/
	/*     HashTable     */
	/*********************/
	
	template<class _Tp>
	inline size_t HashTable<_Tp>::hash(const char* key) const
	{
		size_t h = FIRSTH;
		for(; *key; key++)
			h = (h * A)^(key[0] * B);
		return h % size;
	}
	
	template<class _Tp>
	void HashTable<_Tp>::printHistogram() const
	{
		cout << endl << "Hash Table Contains ";
		cout << getLength() << " Items total" << endl;
		for (int i = 0; i < size; i++)
		{
			cout << i + 1 << ":\t";
			for(int j = 0; j < array[i].getLength(); j++)
				cout << " X";
			cout << endl;
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
	void HashTable<_Tp>::write(ofstream& outF, int fxnNum)
	{
		if(!outF)
			throw runtime_error("Error writing file. Bad ofstream!");
		
		for(size_t i = 0; i < size; i++)
			array[i].write(outF, fxnNum);
	}
}

