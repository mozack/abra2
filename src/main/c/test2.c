#include <iostream>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>

using google::sparse_hash_map;      // namespace where class lives by default
using google::sparse_hash_set;      // namespace where class lives by default
using std::cout;
using std::endl;
using std::tr1::hash;  // or __gnu_cxx::hash, or maybe tr1::hash, depending on your OS

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
	  printf("Comparing: %s : %s\n", s1, s2);
    return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
  }
};

struct my_hash
{
	long operator()(const char* s1) const
	{
		long hash = 0;
		int c;

		while((c = *s1++))
		{
			/* hash = hash * 33 ^ c */
			hash = ((hash << 5) + hash) ^ c;
		}

		return hash;
	}
};

//sparse_hash_map<const char*, int, hash<const char*>, eqstr> months;
sparse_hash_map<const char*, int, my_hash, eqstr> months;

int main()
{
//  sparse_hash_map<const char*, int, hash<const char*>, eqstr> months;

//	cout << "september -> " << months["september"] << endl;
//
//  months["january"] = 31;
//  months["february"] = 28;
//  months["march"] = 31;
//  months["april"] = 30;
//  months["may"] = 31;
//  months["june"] = 30;
//  months["july"] = 31;
//  months["august"] = 31;
//  months["september"] = 30;
//  months["october"] = 31;
//  months["november"] = 30;
//  months["december"] = 31;
//
//  cout << "september -> " << months["september"] << endl;
//  cout << "april     -> " << months["april"] << endl;
//  cout << "june      -> " << months["june"] << endl;
//  cout << "november  -> " << months["november"] << endl;
//
//  char* foo = (char*) malloc(50);
//  char* foo2 = (char*) malloc(50);
//  memset(foo, 0, 50);
//  memset(foo2, 0, 50);
////  memcpy(foo, "september", 9);
//
//  strcpy(foo, "september");
//  strcpy(foo2, foo);
//
//  my_hash H;
//
//  cout << "const hash: " << H("september") << endl;
//  cout << "foo hash: " << H(foo) << endl;
//  cout << "foo2 hash: " << H(foo2) << endl;
//
//  cout << "foo txt: [" << foo << "]" << endl;
//  cout << "foo2 txt: [" << foo2 << "]" << endl;
//
//  cout << "map: " << months[foo] << endl;
//  cout << "map2: " << months["september"] << endl;

	sparse_hash_set<const char*, my_hash, eqstr> Set;

	Set.insert("1234");
}
