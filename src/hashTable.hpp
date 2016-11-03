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

namespace hashTable{
	
	using namespace std;
	
	/******************************/
	/* globally scoped constants */
	/*****************************/
	
	int const HASH_TABLE_SIZE = 20000;
	int const A = 54059;
	int const B = 76963;
	int const C = 86969;
	int const FIRSTH = 37;
	
	/**********************/
	/* class definitions */
	/*********************/
	
	template <class T> class Node;
	template <class T> class LinkedList;
	template <class T> class HashTable;
	
	template <class T>
	class Node {
		friend class LinkedList<T>;
	private:
		Node* next;
		
	public:
		T val;
		
		Node();
	};
	
	template<class T>
	class LinkedList {
	private:
		typename Node<T>::Node* head;
		int length;
		
		void destroyList(typename Node<T>::Node*);
		void insert(const T&, typename Node<T>::Node*);
		typename Node<T>::Node* consolidate(const T&, typename Node<T>::Node*);
		typename Node<T>::Node* getItem(string, typename Node<T>::Node*) const;
		void write(ofstream&, typename Node<T>::Node*);
		void apply(void (T::*fun)(), typename Node<T>::Node*);
		void apply(void (T::*fun)(void*), void*, typename Node<T>::Node*);
		
	public:
		//constructor
		LinkedList();
		~LinkedList();
		
		//modifer
		void insert(const T&);
		typename Node<T>::Node* consolidate(const T&);
		//bool removeItem(string);
		void destroyList();
		void apply(void (T::*fun)());
		void apply(void (T::*fun)(void*), void*);
		
		//properties
		int getLength();
		typename Node<T>::Node* getItem(string) const;
		void write(ofstream&);
		bool empty();
	};
	
	template <class T>
	class HashTable{
	private:
		int size;
		//hash<string> str_hash;
		typename LinkedList<T>::LinkedList* array;
		
		void destroyTable(typename LinkedList<T>::LinkedList*);
		
		size_t hash(const char*) const;
		//size_t hash(string) const;
		
	public:
		//constructor
		HashTable(int s = HASH_TABLE_SIZE);
		~HashTable();
		
		//modifers
		//bool remove(string);
		void destroyTable();
		void insert(const T&, string);
		typename Node<T>::Node* consolidate(const T&, string);
		void apply(void (T::*fun)());
		void apply(void (T::*fun)(void*), void*);
		
		//properties
		int getLength();
		void printTable();
		typename Node<T>::Node* getItem(string) const;
		int getNumItems() const;
		void printHistogram() const;
		void write(ofstream&);
	};
}

#endif
