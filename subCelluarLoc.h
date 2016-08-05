
//#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<cassert>

using namespace std;

struct DBProtein{
	string ID, gene, description, loc;
	
	//constructor
	DBProtein(string);
	DBProtein();
	
	//modifers
	void operator = (const DBProtein&);
	void clear();
	
	//properties
	bool operator < (const DBProtein&) const;
	bool operator > (const DBProtein&) const;
	bool operator == (const DBProtein&) const;
};

DBProtein::DBProtein()
{
	ID = "";
	gene = "";
	description = "";
	loc = "";
}

void DBProtein::clear()
{
	ID = "";
	gene = "";
	description = "";
	loc = "";
	
}

DBProtein::DBProtein(string line)
{
	if(line.empty())
		clear();
	else
	{
		vector<string>elems;
		split(line, '\t', elems);
		assert(elems.size() == 4);
		
		ID = elems[0];
		gene = elems[1];
		description = elems[2];
		loc = elems[3];
	}
}

bool DBProtein::operator < (const DBProtein& compp) const
{
	int comp = strcmp(compp.ID.c_str(), ID.c_str());
	return comp > 0;
}

bool DBProtein::operator > (const DBProtein& compp) const
{
	int comp = strcmp(compp.ID.c_str(), ID.c_str());
	return comp < 0;
}

bool DBProtein::operator == (const DBProtein& compp) const
{
	int comp = strcmp(compp.ID.c_str(), ID.c_str());
	return comp == 0;
}

void DBProtein::operator = (const DBProtein& pget)
{
	ID = pget.ID;
	gene = pget.gene;
	description = pget.description;
	loc = pget.loc;
}

struct Node{
	DBProtein protein;
	Node();
	
	Node *left;
	Node *right;
};

Node::Node()
{
	new struct DBProtein;
	left = nullptr;
	right = nullptr;
}

class Btree{
public:
	Btree();
	~Btree();
	
	void insert(const DBProtein&);
	bool readInProteins(string);
	string locSearch(const DBProtein&) const;
	
private:
	void destroyTree();
	void insert(const DBProtein& p, Node *leaf);
	Node *search(const DBProtein& p, Node *leaf) const;
	Node *search(const DBProtein&) const;
	void destroyTree(Node *leaf);
	
	Node *root;
};

Btree::Btree()
{
	root = nullptr;
}

Btree::~Btree()
{
	destroyTree();
}

void Btree::destroyTree(Node *leaf)
{
	if(leaf != nullptr)
	{
		destroyTree(leaf->left);
		destroyTree(leaf->right);
		delete leaf;
	}
}

void Btree::destroyTree()
{
	destroyTree(root);
}

void Btree::insert(const DBProtein& p, Node *leaf)
{
	if(p < leaf->protein)
	{
		if(leaf->left != nullptr)
			insert(p, leaf->left);
		else
		{
			leaf->left = new Node;
			leaf->left->protein = p;
			leaf->left->left=nullptr;
			leaf->left->right=nullptr;
		}
	}
	else if(p > leaf->protein || p == leaf->protein)
	{
		if(leaf->right != nullptr)
			insert(p, leaf->right);
		else
		{
			leaf->right = new Node;
			leaf->right->protein = p;
			leaf->right->left=nullptr;
			leaf->right->right=nullptr;
		}
	}
}

Node *Btree::search(const DBProtein& p, Node *leaf) const
{
	if(leaf!=nullptr)
	{
		if(p == leaf->protein)
			return leaf;
		if(p < leaf->protein)
			return search(p, leaf->left);
		else
			return search(p, leaf->right);
	}
	else return nullptr;
}

void Btree::insert(const DBProtein& p)
{
	if(root!=nullptr)
		insert(p, root);
	else
	{
		root=new Node;
		root->protein = p;
		root->left=nullptr;
		root->right=nullptr;
	}
}

Node *Btree::search(const DBProtein& p) const
{
	return search(p, root);
}

string Btree::locSearch(const DBProtein& p) const
{
	Node* node = search(p, root);
	if(node == nullptr)
		return "NOT FOUND IN DB";
	
	string loc = node->protein.loc;
	return loc;
}

bool Btree::readInProteins(string fname)
{
	ifstream inF (fname);
	
	if (!inF)
		return false;
	
	string line;
	
	while(!inF.eof()){
		getline(inF, line);
		DBProtein newDBProtein(line);
		if (newDBProtein.ID != "ID") //skip line if it is header line
			insert(newDBProtein);
	}
	
	return true;
}
