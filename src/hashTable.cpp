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
	LinkedList<_Tp>::LinkedList()
	{
		front = nullptr;
		length = 0;
	}
	
	template<class _Tp>
	LinkedList<_Tp>::~LinkedList()
	{
		destroyList();
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::destroyList()
	{
		destroyList(front);
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::destroyList(typename Node<_Tp>::Node* node)
	{
		if(node != nullptr)
		{
			destroyList(node->next);
			delete node;
		}
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::pop_back(const _Tp& newItem)
	{
		/*if(front == nullptr)
		{
			front = new Node<_Tp>;
			front->val = newItem;
			back = front;
			back->prev = front;
			front->next = back;
			front->prev = nullptr;
		}
		else{
			typename Node<_Tp>::Node* temp = new Node<_Tp>;
			temp->val = newItem;
			temp->prev = back;
			temp->next = nullptr;
			back = temp;
			back->prev->next = back;
		}*/
		
		if(front != nullptr)
			pop_back(newItem, front);
		else
		{
			front = new Node<_Tp>;
			front->val = newItem;
			front->next = nullptr;
			length++;
		}
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::pop_back(const _Tp& newItem, typename Node<_Tp>::Node* leaf)
	{
		if(leaf->next != nullptr)
			pop_back(newItem, leaf->next);
		else
		{
			leaf->next = new Node<_Tp>;
			leaf->next->val = newItem;
			leaf->next->next = nullptr;
			length++;
		}
	}
	
	template<class _Tp>
	typename Node<_Tp>::Node* LinkedList<_Tp>::consolidate(const _Tp& newItem)
	{
		if(front != nullptr)
			return consolidate(newItem, front);
		else
		{
			front = new Node<_Tp>;
			front->val = newItem;
			front->next = nullptr;
			length++;
			return front;
		}
	}
	
	template<class _Tp>
	typename Node<_Tp>::Node* LinkedList<_Tp>::consolidate(const _Tp& newItem, typename Node<_Tp>::Node* leaf)
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
			leaf->next = new Node<_Tp>;
			leaf->next->val = newItem;
			leaf->next->next = nullptr;
			length++;
			return leaf->next;
		}
	}

	/*bool LinkedList::removeItem(string key)
	 {
	 if (front->next == nullptr)
		return false;
	 Item * p = front;
	 Item * q = front;
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
	
	template<class _Tp>
	typename Node<_Tp>::Node* LinkedList<_Tp>::getItem(string key, typename Node<_Tp>::Node* node) const
	{
		if(node != nullptr)
		{
			if(node->val == key)
				return node;
			else return getItem(key, node->next);
		}
		else return nullptr;
	}
	
	template<class _Tp>
	typename Node<_Tp>::Node* LinkedList<_Tp>::getItem(string key) const
	{
		return getItem(key, front);
	}
	
	template<class _Tp>
	int LinkedList<_Tp>::getLength()
	{
		return length;
	}
	
	template<class _Tp>
	bool LinkedList<_Tp>::empty()
	{
		return front == nullptr;
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::write(ofstream& outF)
	{
		if(!outF)
			throw runtime_error("Error writing file. Bad ofstream!");
		
		if(front == nullptr)
			return;
		
		write(outF, front);
	}
	
	template<class _Tp>
	void LinkedList<_Tp>::write(ofstream& outF, typename Node<_Tp>::Node* leaf)
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
	
	template<class _Tp>
	void HashTable<_Tp>::destroyTable()
	{
		delete [] array;
	}
	
	template<class _Tp>
	HashTable<_Tp>::~HashTable()
	{
		destroyTable();
	}
	
	template<class _Tp>
	size_t HashTable<_Tp>::hash(const char* key) const
	{
		size_t h = FIRSTH;
		for(; *key; key++)
			h = (h * A) ^ (key[0] * B);
		return h % size;
	}
	
	/*template<class _Tp>
	size_t HashTable<_Tp>::hash(string key) const
	{
		return str_hash(key) % size;
	}*/
	
	template<class _Tp>
	void HashTable<_Tp>::insert(const _Tp& newItem, string key)
	{
		size_t index = hash(key.c_str());
		array[index].pop_back(newItem);
	}
	
	template<class _Tp>
	typename Node<_Tp>::Node* HashTable<_Tp>::consolidate(const _Tp& newItem, string key)
	{
		size_t index = hash(key.c_str());
		return array[index].consolidate(newItem);
	}
	
	/*bool HashTable::remove(string key)
	 {
	 int index = hash(key);
	 return array[index].removeItem(key);
	 }*/
	
	template<class _Tp>
	typename Node<_Tp>::Node* HashTable<_Tp>::getItem(string key) const
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

