//
//  hashTable.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/9/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef hashTable_hpp
#define hashTable_hpp

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <list>

namespace hashTable{
	
	using namespace std;
	
	/******************************/
	/* namespace scoped constants */
	/******************************/
	
	int const DEFAULT_HASH_TABLE_SIZE = 1000;
	int const A = 54059;
	int const B = 76963;
	int const C = 86969;
	int const FIRSTH = 37;
	
	/**********************/
	/* class definitions */
	/*********************/
	
	template<class _Tp> class LinkedList;
	template<class _Tp> class HashTable;
	
	template<class _Tp>
	class LinkedList {
	private:
		list<_Tp> dat;
	
	public:
		//constructor
		LinkedList() {}
		~LinkedList() {}
		
		//modifer
		inline void push_back(const _Tp& newItem){
			dat.push_back(newItem);
		}
		inline void push_front(const _Tp& newItem){
			dat.push_front(newItem);
		}
		_Tp* consolidate(const _Tp&);
		
		//properties
		_Tp* getItem(string);
		void write(ofstream&);
		
		inline size_t getLength() const{
			return dat.size();
		}
		inline bool empty() const{
			return dat.empty();
		}
	};
	
	template<class _Tp>
	class HashTable{
	private:
		const size_t size;
		typename LinkedList<_Tp>::LinkedList* array;
		
		void destroyTable(){
			delete [] array;
		}
		
		inline size_t hash(const char*) const;
		
	public:
		//constructor
		HashTable(size_t s = DEFAULT_HASH_TABLE_SIZE) : size(s){
			array = new LinkedList<_Tp> [size];
		}
		~HashTable(){
			destroyTable();
		}
		
		//modifers
		void insert(const _Tp&, string);
		_Tp* consolidate(const _Tp&, string);
		
		//properties
		int getLength();
		_Tp* getItem(string) const;
		int getNumItems() const;
		void printHistogram() const;
		void write(ofstream&);
	};
}

#endif
