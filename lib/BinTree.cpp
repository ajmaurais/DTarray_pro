//
//  BinTree.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/9/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "BinTree.hpp"

namespace binTree{
	
	template <class _Tp>
	Node<_Tp>::Node()
	{
		left = nullptr;
		right = nullptr;
	}
	
	template <class _Tp>
	_Tp* Node<_Tp>::getVal() const
	{
		return &val;
	}
	
	template<class _Tp> BinTree<_Tp>::BinTree()
	{
		root = nullptr;
	}
	
	template<class _Tp>
	BinTree<_Tp>::~BinTree()
	{
		destroyTree();
	}
	
	template <class _Tp>
	void BinTree<_Tp>::destroyTree(typename Node<_Tp>::Node *leaf)
	{
		if(leaf != nullptr)
		{
			destroyTree(leaf->left);
			destroyTree(leaf->right);
			delete leaf;
		}
	}
	
	template<class _Tp>
	void BinTree<_Tp>::destroyTree()
	{
		destroyTree(root);
	}
	
	template<class _Tp>
	void BinTree<_Tp>::insert(const _Tp& arg, typename Node<_Tp>::Node* leaf)
	{
		if(arg < leaf->val)
		{
			if(leaf->left != nullptr)
				insert(arg, leaf->left);
			else
			{
				leaf->left = new Node<_Tp>;
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
				leaf->right = new Node<_Tp>;
				leaf->right->val = arg;
				leaf->right->left=nullptr;
				leaf->right->right=nullptr;
			}
		}
	}
	
	template <class _Tp>
	typename Node<_Tp>::Node* BinTree<_Tp>::search(const _Tp& arg, typename Node<_Tp>::Node* leaf) const
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
	
	template<class _Tp>
	void BinTree<_Tp>::insert(const _Tp& arg)
	{
		if(root!=nullptr)
			insert(arg, root);
		else
		{
			root = new Node<_Tp>;
			root->val = arg;
			root->left = nullptr;
			root->right = nullptr;
		}
	}
	
	template<class _Tp>
	typename Node<_Tp>::Node* BinTree<_Tp>::search(const _Tp& arg) const
	{
		return search(arg, root);
	}
}
