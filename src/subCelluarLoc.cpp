
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

bool DBProtein::operator == (string comp) const
{
	return comp == ID;
}

void DBProtein::operator = (const DBProtein& pget)
{
	ID = pget.ID;
	gene = pget.gene;
	description = pget.description;
	loc = pget.loc;
}
