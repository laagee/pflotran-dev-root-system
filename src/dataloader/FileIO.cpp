#include "FileIO.h"

FileIO::FileIO(char *filename) {

  printf("%s\n",filename);
  file.open(filename,fstream::in);
  if (!file.is_open()) {
    cout << "ERROR opening file " << filename << ".\n";
    exit(1);
  }
  buffer = NULL;

}

int FileIO::getLine() {

  int ierr = 0;
  delete buffer;

  while(1) {
    string s;
    getline(file,s,'\n');
    if (s.compare(0,1,":") && s.compare(0,1,"!")) {
      buffer = new stringstream(s.c_str());
      break;
    }
  }
  return file.eof() ? 0 : 1;

}

int FileIO::getInputLine() {

  return getLine();

}

int FileIO::readDouble(double *d) {

  *buffer >> *d;
  return buffer->fail() ? 0 : 1;
  
}

int FileIO::readInt(int *i) {

  *buffer >> *i;

  return buffer->fail() ? 0 : 1;
  
}

int FileIO::readWord(char *word) {

  /* Remove any preceding spaces(32), tabs(9), or commas(44) etc */
  char c;
  *buffer >> noskipws >> c;
  while ((c == 32 || c == 44 || c == 9) && c != '\0')
    *buffer >> noskipws >> c;

  /* Copy next group of chars to string */
  string str;
  while (c != 32 && c != 44 && c != 9 && c != '\0') {
    str.append(1,c);
    *buffer >> noskipws >> c;
  }
  strcpy(word,str.c_str());

  /* Remove any trailing spaces or commas etc */
  while ((c == 32 || c == 44 || c == 9) && c != '\0')
    *buffer >> noskipws >> c;
  buffer->unget();

  if (strlen(word) == 0) return 1;
  else return 0;

}

int FileIO::readQuotedWords(char *words) {

  /* Remove any preceding spaces(32), tabs(9), or commas(44) etc */
  char c;
  *buffer >> noskipws >> c;
  while ((c == 32 || c == 44 || c == 9) && c != '\0')
    *buffer >> noskipws >> c;

  string str;
  if (c == 34) { // quote found
    while (c != 34) {
      str.append(1,c);
      *buffer >> noskipws >> c;
    }
  }
  else {
    while (c != 32 && c != 44 && c != 9 && c != '\0') {
      str.append(1,c);
      *buffer >> noskipws >> c;
    }
  }
  strcpy(words,str.c_str());

  /* Remove any trailing spaces or commas etc */
  while ((c == 32 || c == 44 || c == 9) && c != '\0')
    *buffer >> noskipws >> c;
  buffer->unget();

  if (strlen(words) == 0) return 1;
  else return 0;

}

int FileIO::removeQuotes(char *str) {

  /* Remove all quotes */
  string str2;
  while (1) {
    size_t found = str2.find("\"");
    if (found == string::npos) break;
    else str2.erase(found);
  }
  strcpy(str,str2.c_str());

  return 0;

}

int FileIO::findStringInFile(char *card) {

  int ierr = 0;

  file.seekg(0,ios::beg);
  size_t len = strlen(card);
  size_t found = 0;
  while((ierr = getLine()) != 1) {
    string str = buffer->str();
    found = str.find(card);
    if (found!=string::npos){
      break;
    }
  }
  //return found ? 0 : 1;
  return found==string::npos ? 0 : 1;
}

int FileIO::comparesTo(char *str) {
  string str2 = buffer->str();
  return str2.compare(str);
}

int FileIO::startsWith(char *str) {
  return comparesTo(str);
}

void FileIO::checkDefaultMessage(char *word, int *ierr) {

  if (ierr) cout << "\"" << word << "\" set to default value" << endl;
  *ierr = 0;

}

void FileIO::checkErrorMessage(char *word1, char *word2, int ierr) {

  if (ierr) {
    cout << "Error reading \"" << word1 << "\" under keyword \"" << word2 << 
            "\"." << endl;
    exit(1);
  }

}

void FileIO::checkLineErrorMessage(char *word, int ierr) {

  if (ierr) {
    cout << "Error reading in string in \"" << word << "\"." << endl;
    exit(1);
  }

}

void FileIO::toUpper(char *str) {

  int len = (int)strlen(str);
  for (int i=0; i < len; i++) {
    str[i] = toupper(str[i]);
  }

}

void FileIO::toLower(char *str) {

  int len = (int)strlen(str);
  for (int i=0; i < len; i++) {
    str[i] = tolower(str[i]);
  }

}

FileIO::~FileIO() {
  file.close();
}
