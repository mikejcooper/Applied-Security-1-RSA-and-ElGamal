#include "modmul.h"

void filecmp(FILE *fp1, FILE *fp2)
{
	char line1[1000];
    char line2[1000];
    int i = 1;
	while (fgets(line1, sizeof(line1), fp1)) {
		fgets(line2, sizeof(line2), fp2);

		if (strcmp(line1,line2) == 0)
			printf("Test %d: Passed \n", i); 
		else
			printf("Test %d: Failed \n", i); 
		i++;
    }
}

void test(char* file1, char* file2){
	FILE *fp1, *fp2; 
 	if ((fp1 = fopen(file1, "r")) == NULL) {
        printf("cat: can't open %s\n", file1);
    }
    if ((fp2 = fopen(file2, "r")) == NULL) {
        printf("cat: can't open %s\n", file2);
    }
    filecmp(fp1, fp2);
    fclose(fp1);
    fclose(fp2);
}


