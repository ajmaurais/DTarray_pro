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
		typename std::list<_Tp>::iterator it = dat.begin();
		for(;it != dat.end(); it++)
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
		typename list<_Tp>::iterator it = dat.begin();
		for(;it != dat.end(); it++)
			if(*it == key)
				return &(*it);

		return nullptr;
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::write(ofstream& outF)
	{
		if(!outF)
			throw runtime_error("Error writing file. Bad ofstream!");
		
		if(dat.empty())
			return;
		
		typename list<_Tp>::iterator it = dat.begin();
		for(;it != dat.end(); it++)
			(*it).write(outF);
			
	}
	
	/**********************/
	/*     HashTable     */
	/*********************/
	
	template<class _Tp>
	inline size_t HashTable<_Tp>::hash(const char* key) const
	{
		size_t h = FIRSTH;
		for(; *key; key++)
			h = (h * A) ^ (key[0] * B);
		return h % size;
	}
	
	template<class _Tp>
	void HashTable<_Tp>::insert(const _Tp& newItem, string key)
	{
		size_t index = hash(key.c_str());
		array[index].push_front(newItem);
	}
	
	template<class _Tp>
	_Tp* HashTable<_Tp>::consolidate(const _Tp& newItem, string key)
	{
		size_t index = hash(key.c_str());
		return array[index].consolidate(newItem);
	}
	
	template<class _Tp>
	_Tp* HashTable<_Tp>::getItem(string key) const
	{
		size_t index = hash(key.c_str());
		return array[index].getItem(key);
	}
	
	template<class _Tp>
	void HashTable<_Tp>::printHistogram() const
	{
		cout << endl << "Hash Table Contains ";
		cout << getNumItems() << " Items total" << endl;
		for (int i = 0; i < size; i++)
		{
			cout << i + 1 << ":\t";
			for ( int j = 0; j < array[i].getLength(); j++ )
				cout << " X";
			cout << endl;
		}
	}
	
	template<class _Tp>
	int HashTable<_Tp>::getNumItems() const
	{
		int c = 0;
		for (int i = 0; i < size; i++)
			c += array[i].getLength();
		
		return c;
	}
	
	template<class _Tp>
	void HashTable<_Tp>::write(ofstream& outF)
	{
		if(!outF)
			throw runtime_error("Error writing file. Bad ofstream!");
		
		for(int i = 0; i < size; i++)
			array[i].write(outF);
	}
}

