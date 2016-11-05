//
//  linkedList.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/7/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "hashTable.hpp"

namespace hashTable{
	
	/**********************/
	/*        Node        */
	/*********************/
	
	template <class T>
	Node<T>::Node()
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
	void LinkedList<T>::destroyList(typename Node<T>::Node* node)
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
			head = new Node<T>;
			head->val = newItem;
			head->next = nullptr;
			length++;
		}
	}
	
	template<class T>
	void LinkedList<T>::insert(const T& newItem, typename Node<T>::Node* leaf)
	{
		if(leaf->next != nullptr)
			insert(newItem, leaf->next);
		else
		{
			leaf->next = new Node<T>;
			leaf->next->val = newItem;
			leaf->next->next = nullptr;
			length++;
		}
	}
	
	template<class T>
	typename Node<T>::Node* LinkedList<T>::consolidate(const T& newItem)
	{
		if(head != nullptr)
			return consolidate(newItem, head);
		else
		{
			head = new Node<T>;
			head->val = newItem;
			head->next = nullptr;
			length++;
			return head;
		}
	}
	
	template<class T>
	typename Node<T>::Node* LinkedList<T>::consolidate(const T& newItem, typename Node<T>::Node* leaf)
	{
		if(leaf->val == newItem)
		{
			leaf->val.consolidate(newItem);
			return leaf;
		}
		else if(leaf->next != nullptr)
		{
			return consolidate(newItem, leaf->next);
		}
		else
		{
			leaf->next = new Node<T>;
			leaf->next->val = newItem;
			leaf->next->next = nullptr;
			length++;
			return leaf->next;
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
	typename Node<T>::Node* LinkedList<T>::getItem(string key, typename Node<T>::Node* node) const
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
	typename Node<T>::Node* LinkedList<T>::getItem(string key) const
	{
		return getItem(key, head);
	}
	
	template<class T>
	int LinkedList<T>::getLength()
	{
		return length;
	}
	
	template<class T>
	bool LinkedList<T>::empty()
	{
		return head == nullptr;
	}
	
	template<class T>
	void LinkedList<T>::write(ofstream& outF)
	{
		if(!outF)
			throw runtime_error("Error writing file. Bad ofstream!");
		
		if(head == nullptr)
			return;
		
		write(outF, head);
	}
	
	template<class T>
	void LinkedList<T>::write(ofstream& outF, typename Node<T>::Node* leaf)
	{
		if(!outF)
			throw runtime_error("Error writing file. Bad ofstream!");
		
		if(leaf == nullptr)
			return;
		else
		{
			leaf->val.write(outF);
			write(outF, leaf->next);
		}
	}
	
	/**********************/
	/*     HashTable     */
	/*********************/
	
	template<class T>
	HashTable<T>::HashTable(int s)
	{
		size = s;
		array = new LinkedList<T>[DEFAULT_HASH_TABLE_SIZE];
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
	size_t HashTable<T>::hash(const char* key) const
	{
		size_t h = FIRSTH;
		for(; *key; key++)
			h = (h * A) ^ (key[0] * B);
		return h % size;
	}
	
	/*template<class T>
	size_t HashTable<T>::hash(string key) const
	{
		return str_hash(key) % size;
	}*/
	
	template<class T>
	void HashTable<T>::insert(const T& newItem, string key)
	{
		size_t index = hash(key.c_str());
		array[index].insert(newItem);
	}
	
	template<class T>
	typename Node<T>::Node* HashTable<T>::consolidate(const T& newItem, string key)
	{
		size_t index = hash(key.c_str());
		return array[index].consolidate(newItem);
	}
	
	/*bool HashTable::remove(string key)
	 {
	 int index = hash(key);
	 return array[index].removeItem(key);
	 }*/
	
	template<class T>
	typename Node<T>::Node* HashTable<T>::getItem(string key) const
	{
		size_t index = hash(key.c_str());
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
	
	template<class T>
	void HashTable<T>::write(ofstream& outF)
	{
		if(!outF)
			throw runtime_error("Error writing file. Bad ofstream!");
		
		for(int i = 0; i < size; i++)
			array[i].write(outF);
	}
}

