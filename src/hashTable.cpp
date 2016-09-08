//
//  linkedList.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/7/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

/**********************/
/*        Item        */
/*********************/

Item::Item()
{
	Peptide temp;
	val = temp;
}

/**********************/
/*     LinkedList     */
/*********************/

LinkedList::LinkedList()
{
	head = nullptr;
	length = 0;
}

LinkedList::~LinkedList()
{
	destroyList();
}

void LinkedList::destroyList()
{
	destroyList(head);
}

void LinkedList::destroyList(Item* item)
{
	if(item != nullptr)
	{
		destroyList(item->next);
		delete item;
	}
}

void LinkedList::insert(const Peptide& newPeptide)
{
	if(head != nullptr)
		insert(newPeptide, head);
	else
	{
		head = new Item;
		head->val = newPeptide;
		head->next = nullptr;
		length++;
	}
}

void LinkedList::insert(const Peptide& newPeptide, Item* item)
{
	if(item->next != nullptr)
		insert(newPeptide, item->next);
	else
	{
		item->next = new Item;
		item->next->val = newPeptide;
		item->next->next = nullptr;
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

Item* LinkedList::getItem(string key, Item* item) const
{
	if(item != nullptr)
	{
		if(item->val == key)
			return item;
		else return getItem(key, item->next);
	}
	else return nullptr;
}

Item* LinkedList::getItem(string key) const
{
	return getItem(key, head);
}

int LinkedList::getLength()
{
	return length;
}

/**********************/
/*     HashTable     */
/*********************/

HashTable::HashTable(int s)
{
	size = s;
	array = new LinkedList[HASH_TABLE_SIZE];
}

void HashTable::destroyTable()
{
	delete [] array;
}

HashTable::~HashTable()
{
	destroyTable();
}

int HashTable::hash(string key) const
{
	return str_hash(key) % size;
}

void HashTable::insert(const Peptide& newPeptide)
{
	int index = hash(newPeptide.getID());
	array[index].insert(newPeptide);
}

/*bool HashTable::remove(string key)
{
	int index = hash(key);
	return array[index].removeItem(key);
}*/

Item* HashTable::getItem(string itemKey) const
{
	int index = hash(itemKey);
	return array[index].getItem(itemKey);
}

void HashTable::printHistogram() const
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

int HashTable::getNumItems() const
{
	int itemCount = 0;
	for ( int i = 0; i < size; i++ )
		itemCount += array[i].getLength();
	
	return itemCount;
}

