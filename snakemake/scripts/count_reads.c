#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int countReadsInFile(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file: %s\n", filename);
        exit(1);
    }

    int count = 0;
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), file) != NULL) {
        if (buffer[0] == '@') {  // Assumes that each read starts with '@'
            count++;
        }
    }

    fclose(file);
    return count;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: ./lecture_counter <files_list.txt>\n");
        return 1;
    }

    const char *listFilename = argv[1];
    FILE *listFile = fopen(listFilename, "r");
    if (listFile == NULL) {
        printf("Error opening file: %s\n", listFilename);
        return 1;
    }

    int totalLectures = 0;
    char filename[256];
    while (fgets(filename, sizeof(filename), listFile) != NULL) {
        // Remove trailing newline character if present
        if (filename[strlen(filename) - 1] == '\n') {
            filename[strlen(filename) - 1] = '\0';
        }

        int lectures = countReadsInFile(filename);
        printf("%s,%d\n", filename,lectures);
        totalLectures += lectures;
    }

    printf("Total Lectures: %d\n", totalLectures);

    fclose(listFile);
    return 0;
}

