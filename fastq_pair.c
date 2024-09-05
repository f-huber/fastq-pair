/*
 * Pair a fastq file using a quick index.
 *
 * We have two files, and each file holds sequence data. Each data element is represented by only four, and exactly four
 * lines:
 *
 * @id1/1    <- the sequence identifier always begins with an @ sign and then the identifier is the string upto the first whitespace. The identifier must be unique per sequence in a file
 * ATGATC.... <- the DNA sequence
 * +        <- a line that begins with the + sign and may also have the identifier (but doesn't need to!)
 * !%$#^@     <- an ascii representation of the "quality" of the sequence. The quality is a number between 0 and 100 and this is the chr of that number +33
 *
 * In the two files, there should be some sequences that are related, and they are denoted either as 1/2,  forward/reverse, or just by having the whole id the same
 *
 * @id1/1  pairs with @id1/2
 *
 * or
 *
 * @id1/f pairs with @id1/r
 *
 * or
 *
 * @id1 in file 1 pairs with @id1 in file 2
 *
 * We read the file and make a hash of the ID without the /1 or /f and the position of that id in the file (using tell)
 * Then we read the second file and check to see if we have a matching sequence. If we do, we print both sequences
 * one to each file, and we set the "printed" boolean to true in our data structure corresponding to that ID.
 *
 * Finally, we read through our data structure and print out any sequences that we have not seen before.
 *
 * Note that to print out the sequences we seek to the position we recorded and print four lines.
 *
 */

#include "is_gzipped.h"
#include "fastq_pair.h"
#include "robstr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>  // Include zlib for gzip handling

// Function to remove any suffix from a predefined list of possible suffixes
char* removeSuffix(const char* str) {
    // Define the possible suffixes directly within the function
    const char* suffixes[] = {".fastq", ".fastq.gz", ".fq", "fq.gz"};
    int num_suffixes = sizeof(suffixes) / sizeof(suffixes[0]);

    size_t str_len = strlen(str);

    // Loop through all possible suffixes
    for (int i = 0; i < num_suffixes; i++) {
        size_t suffix_len = strlen(suffixes[i]);

        // Check if the string ends with the current suffix
        if (str_len >= suffix_len && strcmp(str + str_len - suffix_len, suffixes[i]) == 0) {
            // Allocate memory for the new string (excluding the suffix)
            char* new_str = (char*)malloc((str_len - suffix_len + 1) * sizeof(char));
            if (new_str == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                exit(1);  // Exit if memory allocation fails
            }
            // Copy the part of the original string excluding the suffix
            strncpy(new_str, str, str_len - suffix_len);
            new_str[str_len - suffix_len] = '\0';  // Null-terminate the new string
            return new_str;
        }
    }
    // If no suffix to remove, return a copy of the original string
    char* new_str = strdup(str);
    return new_str;
}

// Function to write to a regular file or gzip file depending on the flag
void writeToFile(bool is_gzip, gzFile gz_file, FILE* reg_file, const char* line) {
    if (is_gzip) {
        gzprintf(gz_file, "%s", line);  // Write to gzip file
    } else {
        fprintf(reg_file, "%s", line);  // Write to regular file
    }
}

// Function to read a line from either a regular or gzip file
char* readFromFile(bool is_gzip, gzFile gz_file, FILE* reg_file, char* line, int max_len) {
    if (is_gzip) {
        return gzgets(gz_file, line, max_len);  // Read from gzip file
    } else {
        return fgets(line, max_len, reg_file);  // Read from regular file
    }
}

// Void function to seek in either a regular or gzip file
void seekInFile(bool is_gzip, gzFile gz_file, FILE* reg_file, long offset, int whence) {
    if (is_gzip) {
        gzseek(gz_file, offset, whence);  // Seek in gzip file
    } else {
        fseek(reg_file, offset, whence);  // Seek in regular file
    }
}

// Function to get the current position in either a regular or gzip file
long int tellInFile(bool is_gzip, gzFile gz_file, FILE* reg_file) {
    if (is_gzip) {
        return gztell(gz_file);  // Get position in gzip file
    } else {
        return ftell(reg_file);  // Get position in regular file
    }
}

int pair_files(char *left_fn, char *right_fn, struct options *opt) {

    int left_duplicates_counter=0;
    int right_duplicates_counter=0;
    int left_paired_counter=0;
    int right_paired_counter=0;
    int left_single_counter=0;
    int right_single_counter=0;

    // Hash table for the first file (left)
    struct idloc **ids_left;
    ids_left = malloc(sizeof(*ids_left) * opt->tablesize);
    if (ids_left == NULL) {
        fprintf(stderr, "We cannot allocate the memory for a table size of the first file %d. Please try a smaller value for -t\n", opt->tablesize);
        exit(-1);
    }
    // Only allocate memory for ids_right if deduplication is used
    struct idloc **ids_right;
    if (opt->deduplicate) {
        // Hash table for the first file (right)
        ids_right = malloc(sizeof(*ids_right) * opt->tablesize);
        if (ids_right == NULL) {
            fprintf(stderr, "We cannot allocate the memory for a table size of the second file %d. Please try a smaller value for -t\n", opt->tablesize);
            exit(-1);
        }
    }

    FILE *lfp, *rfp;
    gzFile lfp_gz, rfp_gz;
    bool is_gzip_left, is_gzip_right;
    bool is_gzip_out = false;

    is_gzip_left = test_gzip(left_fn);
    is_gzip_right = test_gzip(right_fn);
    if (is_gzip_left || is_gzip_right) {
        is_gzip_out = true;
    }

    fprintf(stderr, "First file is gzipped: %s\n", is_gzip_left ? "true" : "false");
    fprintf(stderr, "Second file is gzipped: %s\n", is_gzip_right ? "true" : "false");
    fprintf(stderr, "Output files will be gzipped: %s\n", is_gzip_out ? "true" : "false");

    char *line = malloc(sizeof(char) * MAXLINELEN + 1);

    if (is_gzip_left){
        if ((lfp_gz = gzopen(left_fn, "rb")) == NULL) {
                fprintf(stderr, "Can't open file %s\n", left_fn);
                exit(1);
        }
    } else {
        if ((lfp = fopen(left_fn, "r")) == NULL) {
                fprintf(stderr, "Can't open file %s\n", left_fn);
                exit(1);
        }
    }
    char *aline; /* this variable is not used, it suppresses a compiler warning */

    long int nextposition = 0;

    /*
     * Read the first file and make an index of that file.
     */
    while (1) {
        aline = readFromFile(is_gzip_left, lfp_gz, lfp, line, MAXLINELEN);
        if (aline == NULL) {
            break;  // End of file
        }

        struct idloc *newid;
        newid = (struct idloc *) malloc(sizeof(*newid));

        if (newid == NULL) {
            fprintf(stderr, "Can't allocate memory for new ID pointer - first file\n");
            return 0;
        }

        line[strcspn(line, "\n")] = '\0';
        if (opt->splitspace)
            line[strcspn(line, " \t")] = '\0';

        /*
         * Figure out what the match mechanism is. We have four examples so
         *     i.   using /1 and /2
         *     ii.  using /f and /r
	     *     iii. using ' 1...' and ' 2....'
         *     iii. just having the whole name
         *
         * If there is a /1 or /2 in the file name, we set that part to null so the string is only up
         * to before the / and use that to store the location.
         */

        char lastchar = line[strlen(line)-1];
        char lastbutone = line[strlen(line)-2];
        if ('/' == lastbutone || '_' == lastbutone || '.' == lastbutone){
            if ('1' == lastchar || '2' == lastchar || 'f' == lastchar ||  'r' == lastchar){
                line[strlen(line)-1] = '\0'; // Add the null terminator at the new end of the string
            }
        } else {
            line[strlen(line)+1] = '\0';
            line[strlen(line)-1] = '/';
        }

        if (opt->verbose)
            fprintf(stderr, "ID first file is |%s|\n", line);

        // Hash the ID
        unsigned hashval = hash(line) % opt->tablesize;

        // Check if the ID already exists in the hash table (duplicate)
        if (opt->deduplicate) {
            struct idloc *existing = ids_left[hashval];
            while (existing != NULL) {
                if (strcmp(existing->id, line) == 0) {
                    // ID already exists, do not add it again
                    if (opt->verbose)
                        fprintf(stderr, "Duplicate ID found in the first file, skipping: %s\n", line);
                        left_duplicates_counter++;
                    free(newid);  // Free the memory we allocated for newid
                    newid = NULL;
                    break;
                }
                existing = existing->next;
            }
        }


        if (newid != NULL) {
            // If the ID is not a duplicate, proceed with adding it to the hash table
            newid->id = dupstr(line);
            newid->pos = nextposition;
            newid->printed = false;
            newid->next = ids_left[hashval];  // Insert at the head of the list
            ids_left[hashval] = newid;
        }

        /* read the next three lines and ignore them: sequence, header, and quality */
        for (int i=0; i<3; i++)
            aline = readFromFile(is_gzip_left, lfp_gz, lfp, line, MAXLINELEN);

        // Get the current position using tellInFile
        nextposition = tellInFile(is_gzip_left, lfp_gz, lfp);
    }


    /*
     * Now just print all the id lines and their positions
     */

    if (opt->print_table_counts) {
        fprintf(stdout, "Bucket sizes\n");

        for (int i = 0; i < opt->tablesize; i++) {
            struct idloc *ptr = ids_left[i];
            int counter=0;
            while (ptr != NULL) {
                // fprintf(stdout, "ID: %s Position %ld\n", ptr->id, ptr->pos);
                counter++;
                ptr = ptr->next;
            }
            fprintf(stdout, "%d\t%d\n", i, counter);
        }
    }

   /* now we want to open output files for left_paired, right_paired, and right_single */

    FILE *left_paired, *left_single, *right_paired, *right_single;
    gzFile left_paired_gz, left_single_gz, right_paired_gz, right_single_gz;


    char *lpfn, *rpfn, *lsfn, *rsfn;
    if (is_gzip_out) {
        lpfn = catstr(removeSuffix(left_fn), ".paired.fastq.gz");
        rpfn = catstr(removeSuffix(right_fn), ".paired.fastq.gz");
        lsfn = catstr(removeSuffix(left_fn), ".single.fastq.gz");
        rsfn = catstr(removeSuffix(right_fn), ".single.fastq.gz");
    } else {
        lpfn = catstr(removeSuffix(left_fn), ".paired.fastq");
        rpfn = catstr(removeSuffix(right_fn), ".paired.fastq");
        lsfn = catstr(removeSuffix(left_fn), ".single.fastq");
        rsfn = catstr(removeSuffix(right_fn), ".single.fastq");
    }

    printf("Writing the paired reads to %s and %s\nWriting the single reads to %s and %s\n", lpfn, rpfn, lsfn, rsfn);

    // Create output files
    if (is_gzip_out){
        if ((left_paired_gz = gzopen(lpfn, "wb")) == NULL) {
                fprintf(stderr, "Can't open file %s\n", lpfn);
                exit(1);
        }

        if ((left_single_gz = gzopen(lsfn, "wb")) == NULL) {
            fprintf(stderr, "Can't open file %s\n", lsfn);
            exit(1);
        }

        if ((right_paired_gz = gzopen(rpfn, "wb")) == NULL) {
            fprintf(stderr, "Can't open file %s\n", rpfn);
            exit(1);
        }

        if ((right_single_gz = gzopen(rsfn, "wb")) == NULL) {
            fprintf(stderr, "Can't open file %s\n", rsfn);
            exit(1);
        }    

    } else {
        if ((left_paired = fopen(lpfn, "w")) == NULL ) {
            fprintf(stderr, "Can't open file %s\n", lpfn);
            exit(1);
        }

        if ((left_single = fopen(lsfn, "w")) == NULL) {
            fprintf(stderr, "Can't open file %s\n", lsfn);
            exit(1);
        }

        if ((right_paired = fopen(rpfn, "w")) == NULL) {
            fprintf(stderr, "Can't open file %s\n", rpfn);
            exit(1);
        }

        if ((right_single = fopen(rsfn, "w")) == NULL) {
            fprintf(stderr, "Can't open file %s\n", rsfn);
            exit(1);
        }        
    }
    /*
    * Now read the second file, and print out things in common
    */

    if (is_gzip_right){
        if ((rfp_gz = gzopen(right_fn, "rb")) == NULL) {
                fprintf(stderr, "Can't open file %s\n", left_fn);
                exit(1);
        }
    } else {
        if ((rfp = fopen(right_fn, "r")) == NULL) {
                fprintf(stderr, "Can't open file %s\n", left_fn);
                exit(1);
        }
    }

    nextposition = 0;

    while (1) {
        aline = readFromFile(is_gzip_right, rfp_gz, rfp, line, MAXLINELEN);

        if (aline == NULL) {
            break;  // End of file
        }

        struct idloc *newid;
        newid = (struct idloc *) malloc(sizeof(*newid));

        if (newid == NULL) {
            fprintf(stderr, "Can't allocate memory for new ID pointer - second file\n");
            return 0;
        }

        // make a copy of the current line so we can print it out later.
        char *headerline = dupstr(line);

        line[strcspn(line, "\n")] = '\0';
        if (opt->splitspace)
            line[strcspn(line, " \t")] = '\0';

        /* remove the last character, as we did above */
        char lastchar = line[strlen(line)-1];
        char lastbutone = line[strlen(line)-2];
        if ('/' == lastbutone || '_' == lastbutone || '.' == lastbutone){
            if ('1' == lastchar || '2' == lastchar || 'f' == lastchar ||  'r' == lastchar){
                line[strlen(line)-1] = '\0'; // Add the null terminator at the new end of the string
            }
        } else {
            line[strlen(line)+1] = '\0';
            line[strlen(line)-1] = '/';
        }

        if (opt->verbose)
            fprintf(stderr, "ID second file is |%s|\n", line);

        // Hash the ID
        unsigned hashval = hash(line) % opt->tablesize;

        // Store the current identifier outside of the line variable
        char * entryid = dupstr(line);

        // Check if the ID already exists in the hash table (duplicate)
        if (opt->deduplicate) {
            struct idloc *existing = ids_right[hashval];
            while (existing != NULL) {
                if (strcmp(existing->id, line) == 0) {
                    // ID already exists, do not add it again
                    if (opt->verbose)
                        fprintf(stderr, "Duplicate ID found in the second file, skipping: %s\n", line);
                        right_duplicates_counter++;
                    free(newid);  // Free the memory we allocated for newid
                    newid = NULL;
                    break;
                }
                existing = existing->next;
            }
        }

        if (newid != NULL) {
            if (opt->deduplicate) {
                // If the ID is not a duplicate, proceed with adding it to the hash table of the second file
                newid->id = dupstr(line);
                newid->pos = nextposition;
                newid->printed = false;
                newid->next = ids_right[hashval];  // Insert at the head of the list
                ids_right[hashval] = newid;
            }
            // now see if we have the mate pair
            unsigned hashval = hash(line) % opt->tablesize;
            struct idloc *ptr = ids_left[hashval];
            long int posn = -1; // -1 is not a valid file position
            while (ptr != NULL) {
                if (strcmp(ptr->id, line) == 0) {
                    posn = ptr->pos;
                    ptr->printed = true;
                }
                ptr = ptr->next;
            }

            if (posn != -1) {
                // we have a match.
                // lets process the left file
                seekInFile(is_gzip_left, lfp_gz, lfp, posn, SEEK_SET);
                left_paired_counter++;
                for (int i=0; i<=3; i++) {
                    aline = readFromFile(is_gzip_left, lfp_gz, lfp, line, MAXLINELEN);
                    if (i == 0 && opt->formatid) {
                        writeToFile(is_gzip_out, left_paired_gz, left_paired, catstr(entryid, "1\n"));
                    } else {
                        writeToFile(is_gzip_out, left_paired_gz, left_paired, line);
                    }
                }
                // now process the right file
                if (opt->formatid) {
                    writeToFile(is_gzip_out, right_paired_gz, right_paired, catstr(entryid, "2\n"));
                } else {
                    writeToFile(is_gzip_out, right_paired_gz, right_paired, headerline);
                }
                right_paired_counter++;
                for (int i=0; i<=2; i++) {
                    aline = readFromFile(is_gzip_right, rfp_gz, rfp, line, MAXLINELEN);
                    writeToFile(is_gzip_out, right_paired_gz, right_paired, line);
                }
            }
            else {
                if (opt->formatid) {
                    writeToFile(is_gzip_out, right_single_gz, right_single, catstr(entryid, "2\n"));
                } else {
                    writeToFile(is_gzip_out, right_single_gz, right_single, headerline);
                }
                right_single_counter++;
                for (int i=0; i<=2; i++) {
                    aline = readFromFile(is_gzip_right, rfp_gz, rfp, line, MAXLINELEN);
                    writeToFile(is_gzip_out, right_single_gz, right_single, line);
                }
            }
        } else {
            for (int i=0; i<=2; i++) {
                aline = readFromFile(is_gzip_right, rfp_gz, rfp, line, MAXLINELEN);
            }
        }
    }

    /* all that remains is to print the unprinted singles from the left file */

    for (int i = 0; i < opt->tablesize; i++) {
        struct idloc *ptr = ids_left[i];
        while (ptr != NULL) {
            if (! ptr->printed) {
                seekInFile(is_gzip_left, lfp_gz, lfp, ptr->pos, SEEK_SET);
                left_single_counter++;
                for (int n=0; n<=3; n++) {
                    aline = readFromFile(is_gzip_left, lfp_gz, lfp, line, MAXLINELEN);
                    if (n == 0 && opt->formatid) {
                        writeToFile(is_gzip_out, left_single_gz, left_single, catstr(ptr->id, "1\n"));
                    } else {
                        writeToFile(is_gzip_out, left_single_gz, left_single, line);
                    }
                }
            }
            ptr = ptr->next;
        }
    }

    fprintf(stdout, "Left paired: %-14d Right paired: %d \nLeft single: %-14d Right single: %d\n",
            left_paired_counter, right_paired_counter, left_single_counter, right_single_counter);
    if (opt->deduplicate) {
        fprintf(stdout, "Left duplicates: %-10d Right duplicates: %d\n",
                left_duplicates_counter, right_duplicates_counter);
    }

    if (is_gzip_left) {
        gzclose(lfp_gz);
    } else {
        fclose(lfp);
    }

    if (is_gzip_right) {
        gzclose(rfp_gz);
    } else {
        fclose(rfp);
    }

    if (is_gzip_out) {
        gzclose(left_paired_gz);
        gzclose(left_single_gz);
        gzclose(right_paired_gz);
        gzclose(right_single_gz);
    } else {
        fclose(left_paired);
        fclose(left_single);
        fclose(right_paired);
        fclose(right_single);
    }

    /*
     * Free up the memory for all the pointers
     */


    for (int i = 0; i < opt->tablesize; i++) {
        struct idloc *ptr = ids_left[i];
        struct idloc *next;
        while (ptr != NULL) {
            next = ptr->next;
            free(ptr);
            ptr=next;
        }
    }

    free(ids_left);
    free(line);

    if (opt->deduplicate) {
        for (int i = 0; i < opt->tablesize; i++) {
            struct idloc *ptr = ids_right[i];
            struct idloc *next;
            while (ptr != NULL) {
                next = ptr->next;
                free(ptr);
                ptr=next;
            }
        }
        free(ids_right);
    }

    return 0;
}


unsigned hash (char *s) {
    unsigned hashval;

    for (hashval=0; *s != '\0'; s++)
        hashval = *s + 31 * hashval;
    return hashval;
}
