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
		std::list<_Tp> dat;
	
	public:
		//constructor
		LinkedList() {}
		~LinkedList() {}
		
		//modifer
		void push_back(const _Tp& newItem){
			dat.push_back(newItem);
		}
		void push_front(const _Tp& newItem){
			dat.push_front(newItem);
		}
		_Tp* consolidate(const _Tp&);
		_Tp* getItem(std::string);
		void apply(int fxnNum = 0);
		void write(std::ofstream& outF, int fxnNum = 0);
		
		//properties
		bool itemExists(std::string) const;
		size_t getLength() const{
			return dat.size();
		}
		bool empty() const{
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
		
		inline size_t hash(std::string) const;
		
	public:
		//constructor
		HashTable(size_t s = DEFAULT_HASH_TABLE_SIZE) : size(s){
			array = new LinkedList<_Tp> [size];
		}
		~HashTable(){
			destroyTable();
		}
		
		//modifers
		void insert(const _Tp& newItem, std::string key){
			array[hash(key)].push_front(newItem);
		}
		_Tp* consolidate(const _Tp& newItem, std::string key){
			return array[hash(key)].consolidate(newItem);
		}
		
		//properties
		bool itemExists(std::string key) const{
			return array[hash(key)].itemExists(key);
		}
		_Tp* getItem(std::string key) const {
			return array[hash(key)].getItem(key);
		}
		unsigned long getLength() const;
		void printHistogram() const; //for debuging
		void write(std::ofstream& outF, int fxnNum = 0);
		void apply(int fxnNum = 0);
	};//end of HashTable class
}//end of hashTable namespace

#endif
