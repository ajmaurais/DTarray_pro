//
//  BinTree.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/9/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include<iostream>

using namespace std;

namespace binTree{
	
	/**********************/
	/* class definitions */
	/*********************/
	
	template <class T> class BinTree;
	template <class T> class Node;
	
	template <class T>
	class Node{
		friend class BinTree<T>;
	public:
		T val;
		
		Node();
		
		T getVal() const;
	private:
		Node *left;
		Node *right;
	};
	
	template <class T>
	class BinTree{
	public:
		//constructor
		BinTree();
		~BinTree();
		
		//modifers
		void insert(const T&);
		void destroyTree();
		
		//properties
		typename Node<T>::Node* search(const T&) const;
		
	private:
		typename Node<T>::Node* root;
		
		//modifers
		void insert(const T&, typename Node<T>::Node*);
		void destroyTree(typename Node<T>::Node*);
		
		//properties
		typename Node<T>::Node* search(const T&, typename Node<T>::Node*) const;
	};
}
