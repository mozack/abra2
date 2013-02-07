#include <iostream>
#include <sparsehash/sparse_hash_map>

using google::sparse_hash_map;      // namespace where class lives by default
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

sparse_hash_map<const char*, int, hash<const char*>, eqstr> months;

int main()
{
//  sparse_hash_map<const char*, int, hash<const char*>, eqstr> months;

	cout << "september -> " << months["september"] << endl;

  months["january"] = 31;
  months["february"] = 28;
  months["march"] = 31;
  months["april"] = 30;
  months["may"] = 31;
  months["june"] = 30;
  months["july"] = 31;
  months["august"] = 31;
  months["september"] = 30;
  months["october"] = 31;
  months["november"] = 30;
  months["december"] = 31;

  cout << "september -> " << months["september"] << endl;
  cout << "april     -> " << months["april"] << endl;
  cout << "june      -> " << months["june"] << endl;
  cout << "november  -> " << months["november"] << endl;

  char* foo = (char*) malloc(50);
  char* foo2 = (char*) malloc(50);
  memset(foo, 0, 50);
  memset(foo2, 0, 50);
//  memcpy(foo, "september", 9);

  strcpy(foo, "september");
  strcpy(foo2, foo);

  hash<char*> H;

  cout << "const hash: " << H("september") << endl;
  cout << "foo hash: " << H(foo) << endl;
  cout << "foo2 hash: " << H(foo2) << endl;

  cout << "foo txt: [" << foo << "]" << endl;
  cout << "foo2 txt: [" << foo2 << "]" << endl;

  cout << "map: " << months[foo] << endl;
  cout << "map2: " << months["september"] << endl;
}
