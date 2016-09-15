//
//  linkedList.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/7/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "hashTable.hpp"

/**********************/
/*        Item        */
/*********************/

template <class T>
LNode<T>::LNode()
{
	next = nullptr;
}

/**********************/
/*     LinkedList     */
/*********************/

template<class T>
LinkedList<T>::LinkedList()
{
	head = nullptr;
	length = 0;
}

template<class T>
LinkedList<T>::~LinkedList()
{
	destroyList();
}

template<class T>
void LinkedList<T>::destroyList()
{
	destroyList(head);
}

template<class T>
void LinkedList<T>::destroyList(typename LNode<T>::LNode* node)
{
	if(node != nullptr)
	{
		destroyList(node->next);
		delete node;
	}
}

template<class T>
void LinkedList<T>::insert(const T& newItem)
{
	if(head != nullptr)
		insert(newItem, head);
	else
	{
		head = new LNode<T>;
		head->val = newItem;
		head->next = nullptr;
		length++;
	}
}

template<class T>
void LinkedList<T>::insert(const T& newItem, typename LNode<T>::LNode* leaf)
{
	if(leaf->next != nullptr)
		insert(newItem, leaf->next);
	else
	{
		leaf->next = new LNode<T>;
		leaf->next->val = newItem;
		leaf->next->next = nullptr;
		length++;
	}
}

/*bool LinkedList::removeItem(string key)
{
	if (head->next == nullptr)
		return false;
	Item * p = head;
	Item * q = head;
	while (q != nullptr)
	{
		if (q->val == key)
		{
			p->next = q->next;
			delete q;
			length--;
			return true;
		}
		p = q;
		q = p->next;
	}
	return false;
}*/

template<class T>
typename LNode<T>::LNode* LinkedList<T>::getItem(string key, typename LNode<T>::LNode* node) const
{
	if(node != nullptr)
	{
		if(node->val == key)
			return node;
		else return getItem(key, node->next);
	}
	else return nullptr;
}

template<class T>
typename LNode<T>::LNode* LinkedList<T>::getItem(string key) const
{
	return getItem(key, head);
}

template<class T>
int LinkedList<T>::getLength()
{
	return length;
}

/**********************/
/*     HashTable     */
/*********************/

template<class T>
HashTable<T>::HashTable(int s)
{
	size = s;
	array = new LinkedList<T>[HASH_TABLE_SIZE];
}

template<class T>
void HashTable<T>::destroyTable()
{
	delete [] array;
}

template<class T>
HashTable<T>::~HashTable()
{
	destroyTable();
}

template<class T>
int HashTable<T>::hash(string key) const
{
	return str_hash(key) % size;
}

template<class T>
void HashTable<T>::insert(const T& newItem, string key)
{
	int index = hash(key);
	array[index].insert(newItem);
}

/*bool HashTable::remove(string key)
{
	int index = hash(key);
	return array[index].removeItem(key);
}*/

template<class T>
typename LNode<T>::LNode* HashTable<T>::getItem(string key) const
{
	int index = hash(key);
	return array[index].getItem(key);
}

template<class T>
void HashTable<T>::printHistogram() const
{
	cout << "\nHash Table Contains ";
	cout << getNumItems() << " Items total\n";
	for ( int i = 0; i < size; i++ )
	{
		cout << i + 1 << ":\t";
		for ( int j = 0; j < array[i].getLength(); j++ )
			cout << " X";
		cout << "\n";
	}
}

template<class T>
int HashTable<T>::getNumItems() const
{
	int itemCount = 0;
	for ( int i = 0; i < size; i++ )
		itemCount += array[i].getLength();
	
	return itemCount;
}

