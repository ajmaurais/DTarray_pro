//
//  BinTree.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/9/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "BinTree.hpp"

template <class T>
Node<T>::Node()
{
	//val = new T;
	left = nullptr;
	right = nullptr;
}

template <class T>
T Node<T>::getVal() const
{
	return val;
}

template<class T> BinTree<T>::BinTree()
{
	root = nullptr;
}

template<class T>
BinTree<T>::~BinTree()
{
	destroyTree();
}

template <class T>
void BinTree<T>::destroyTree(typename Node<T>::Node *leaf)
{
	if(leaf != nullptr)
	{
		destroyTree(leaf->left);
		destroyTree(leaf->right);
		delete leaf;
	}
}

template<class T>
void BinTree<T>::destroyTree()
{
	destroyTree(root);
}

template<class T>
void BinTree<T>::insert(const T& arg, typename Node<T>::Node* leaf)
{
	if(arg < leaf->val)
	{
		if(leaf->left != nullptr)
			insert(arg, leaf->left);
		else
		{
			leaf->left = new Node<T>;
			leaf->left->val = arg;
			leaf->left->left=nullptr;
			leaf->left->right=nullptr;
		}
	}
	else if(arg > leaf->val || arg == leaf->val)
	{
		if(leaf->right != nullptr)
			insert(arg, leaf->right);
		else
		{
			leaf->right = new Node<T>;
			leaf->right->val = arg;
			leaf->right->left=nullptr;
			leaf->right->right=nullptr;
		}
	}
}

template <class T>
typename Node<T>::Node* BinTree<T>::search(const T& arg, typename Node<T>::Node* leaf) const
{
	if(leaf!=nullptr)
	{
		if(arg == leaf->val)
			return leaf;
		if(arg < leaf->val)
			return search(arg, leaf->left);
		else
			return search(arg, leaf->right);
	}
	else return nullptr;
}

template<class T>
void BinTree<T>::insert(const T& arg)
{
	if(root!=nullptr)
		insert(arg, root);
	else
	{
		root = new Node<T>;
		root->val = arg;
		root->left = nullptr;
		root->right = nullptr;
	}
}

template<class T>
typename Node<T>::Node* BinTree<T>::search(const T& arg) const
{
	return search(arg, root);
}
