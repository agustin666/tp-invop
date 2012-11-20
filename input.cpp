#include <cstring>
#include <cstdio>
#include <string>
using namespace std;

char* getParam(string needle, char* haystack[], int count) {
	int i = 0;
	for (i = 0; i < count; i ++) {
		if (strcmp(needle.c_str(), haystack[i]) == 0) {
			if (i < count -1) {
				return haystack[i+1];
			}	
		}
	}
	return 0;	
}


int isParam(string needle, char* haystack[], int count) {
	int i = 0;
	for (i = 0; i < count; i ++) {
		if (strcmp(needle.c_str(), haystack[i]) == 0) {
			return 1;
		}
	}
	return 0;	
}

