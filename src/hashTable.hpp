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
	
	template<class _Tp> class Node;
	template<class _Tp> class LinkedList;
	template<class _Tp> class HashTable;
	
	template<class _Tp>
	class Node {
		friend class LinkedList<_Tp>;
	private:
		Node* next;
		Node* prev;
		
	public:
		_Tp val;
		
		Node(){
			next = nullptr;
			prev = nullptr;
		}
	};
	
	template<class _Tp>
	class LinkedList {
	private:
		typename Node<_Tp>::Node* front;
		//typename Node<_Tp>::Node* back;
		int length;
		
		void destroyList(typename Node<_Tp>::Node*);
		void pop_back(const _Tp&, typename Node<_Tp>::Node*);
		//void pop_front(const _Tp&, typename Node<_Tp>::Node*);
		typename Node<_Tp>::Node* consolidate(const _Tp&, typename Node<_Tp>::Node*);
		typename Node<_Tp>::Node* getItem(string, typename Node<_Tp>::Node*) const;
		void write(ofstream&, typename Node<_Tp>::Node*);
		
	public:
		//constructor
		LinkedList();
		~LinkedList();
		
		/*class iterator{
		private:
			typename Node<_Tp>::Node* it;
		public:
			void operator ++ ();
			
		};*/
		
		//modifer
		void pop_back(const _Tp&);
		void pop_front(const _Tp&);
		typename Node<_Tp>::Node* consolidate(const _Tp&);
		//bool removeItem(string);
		void destroyList();
		
		//properties
		int getLength();
		typename Node<_Tp>::Node* getItem(string) const;
		void write(ofstream&);
		bool empty();
	};
	
	template<class _Tp>
	class HashTable{
	private:
		const size_t size;
		//hash<string> str_hash;
		typename LinkedList<_Tp>::LinkedList* array;
		
		void destroyTable(typename LinkedList<_Tp>::LinkedList*);
		
		size_t hash(const char*) const;
		//size_t hash(string) const;
		
	public:
		//constructor
		HashTable(size_t s = DEFAULT_HASH_TABLE_SIZE) : size(s){
			array = new LinkedList<_Tp> [size];
		}
		~HashTable();
		
		//modifers
		//bool remove(string);
		void destroyTable();
		void insert(const _Tp&, string);
		typename Node<_Tp>::Node* consolidate(const _Tp&, string);
		
		//properties
		int getLength();
		void printTable();
		typename Node<_Tp>::Node* getItem(string) const;
		int getNumItems() const;
		void printHistogram() const;
		void write(ofstream&);
	};
}

#endif
