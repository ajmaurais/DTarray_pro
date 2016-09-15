//
//  Header.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/9/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include <functional>
#include <iostream>

using namespace std;

/******************************/
/* globally scoped constants */
/*****************************/

int const HASH_TABLE_SIZE = 20000;

/**********************/
/* class definitions */
/*********************/

template <class T> class LNode;
template <class T> class LinkedList;
template <class T> class HashTable;

template <class T>
class LNode {
	friend class LinkedList<T>;
private:
	LNode* next;

public:
	T val;
	
	LNode();
};

template<class T>
class LinkedList {
private:
	typename LNode<T>::LNode* head;
	int length;
	
	void destroyList(typename LNode<T>::LNode*);
	void insert(const T&, typename LNode<T>::LNode*);
	typename LNode<T>::LNode* getItem(string, typename LNode<T>::LNode*) const;
	
public:
	//constructor
	LinkedList();
	~LinkedList();
	
	//modifer
	void insert(const T&);
	//bool removeItem(string);
	void destroyList();
	void 
	
	//properties
	int getLength();
	typename LNode<T>::LNode* getItem(string) const;
};

template <class T>
class HashTable{
private:
	int size;
	typename LinkedList<T>::LinkedList* array;
	hash<string> str_hash;
	
	int hash(string) const;
	void destroyTable(typename LinkedList<T>::LinkedList*);
public:
	//constructor
	HashTable(int s = HASH_TABLE_SIZE);
	~HashTable();
	
	//modifers
	//bool remove(string);
	void destroyTable();
	void insert(const T&, string);
	
	//properties
	int getLength();
	void printTable();
	typename LNode<T>::LNode* getItem(string) const;
	int getNumItems() const;
	void printHistogram() const;
};
