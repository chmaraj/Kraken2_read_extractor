import timeit
import multiprocessing
import gzip


class ExtractReads(object):

    def __init__(self, args):
        self.args = args
        self.input = args.input
        self.reads = args.reads
        self.selection = args.selection
        self.output = args.output
        self.threads = int(args.threads)
        self.all_subtaxa = args.all_subtaxa
        if args.all_subtaxa:
            if not args.kraken_report:
                print("Kraken2 report must be provided if using --all_subtaxa flag")
                quit()
            else:
                self.kraken_report = args.kraken_report
                self.all_classifications = []
        self.read_ids = []

        self.run()

    def run(self):
        start = timeit.default_timer()
        print("Finding read IDs associated with classification")
        if self.all_subtaxa:
            self.identify_classifications()

            # Split kraken output into [threads] number of equally sized bins of lines for multiprocessing
            with open(self.input, "r") as f:
                lines = f.read().split("\n")
                segmented_lines = []
                for i in range(self.threads):
                    segment = lines[i::self.threads]
                    segmented_lines.append(segment)
                pool = multiprocessing.Pool(self.threads)
                results = pool.map_async(self.identify_subtaxa_ids, segmented_lines)
                pool.close()
                pool.join()

                # Join together results into single list containing all read ids associated with subtaxa classifications
                all_results = results.get()
                for result in all_results:
                    for read_id in result:
                        if read_id not in self.read_ids:
                            self.read_ids.append(read_id)
        else:
            with open(self.input, "r") as f:
                lines = f.read().split("\n")
                for line in lines:
                    if not line:
                        continue
                    fields = line.split("\t")
                    if self.selection in fields[2]:
                        self.read_ids.append(fields[1])

        print("Chunking read file for parallel processing")
        output_file = "{}/{}_selected_reads.fastq".format(self.output, self.selection)
        with open(output_file, "w") as out_file:
            lines = []
            if self.reads.split(".")[-1] == "gz":
                with gzip.open(self.reads, "rt") as f:
                    lines = f.read().strip().split("\n")
            else:
                with open(self.reads, "r") as f:
                    lines = f.read().strip().split("\n")
            all_entries = []

            # Split and parse fastq file. Much faster than using SeqIO.parse
            for i in range(int(len(lines) / 4)):
                entry_start = i * 4
                entry = "\n".join(lines[entry_start:entry_start + 4])
                all_entries.append(entry)
            segmented_entries = []

            # Split reads into [threads] number of equally sized bins for multiprocessing
            for i in range(self.threads):
                segment = all_entries[i::self.threads]
                segmented_entries.append(segment)

            print("Extracting reads")
            pool = multiprocessing.Pool(self.threads)
            results = pool.map_async(self.locate_reads, segmented_entries)
            pool.close()
            pool.join()

            # Combine all located reads into single variable and then write to output file
            all_results = results.get()
            for result in all_results:
                for entry in result:
                    out_file.write("{}\n".format(entry))
            end = timeit.default_timer()
            print("Completed in {} seconds".format(str(end-start)))

    def identify_classifications(self):
        # Identify all subtaxa based on parsing the kraken report.
        # Subtaxa are identified as those classifications listed under the primary classifier which are indented
        # more than the primary classifier in the report, but are before the next classification at the same
        # indentation level as the primary classifier.
        classification_located = False
        classification_index = 0
        group_end_index = 0
        with open(self.kraken_report, "r") as f:
            lines = f.read().strip().split("\n")
            i = 0
            # First step is to find the line with the primary classifier
            while not classification_located:
                linecount = 0
                for line in lines:
                    if not line:
                        continue
                    fields = line.split("\t")
                    if fields[-1][i:] == self.selection:
                        classification_located = True
                        classification_index = linecount
                        break
                    linecount += 1
                i += 1

            # Next find where the next line is with the same indentation level as the primary classifier
            for j in range(classification_index + 1, len(lines)):
                fields = lines[j].split("\t")
                if not fields[-1][i:].startswith(" "):
                    group_end_index = j
                    break
                if j == len(lines):
                    group_end_index = j

            # All classifications between the primary classifier and the end point are considered subtaxa
            for index in range(classification_index, group_end_index):
                fields = lines[index].split("\t")
                group = fields[-1].strip()
                if group not in self.all_classifications:
                    self.all_classifications.append(group)

    def identify_subtaxa_ids(self, lines):
        # Find read ids associated with the subtaxa classifications from the kraken output
        return_ids = []
        for line in lines:
            if not line:
                continue
            fields = line.split("\t")
            read_id = fields[1]
            read_classification = fields[2].split("(")[0].strip()
            for item in self.all_classifications:
                if item == read_classification:
                    if read_id not in return_ids:
                        return_ids.append(read_id)
        return return_ids

    def locate_reads(self, entries_list):
        # Find reads based on list of read ids and return the full list once completed
        located_reads = []
        for entry in entries_list:
            entry_id = entry.split("\n")[0]  # only care about read ID in fastq entry for the sake of locating them
            if any(read in entry_id for read in self.read_ids):
                located_reads.append(entry)
        return located_reads


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Extract reads based on Kraken2 classification')

    parser.add_argument('-i', '--input', required=True,
                        help='The output file from Kraken2')
    parser.add_argument('-r', '--reads', required=True,
                        help='The file to extract reads from')
    parser.add_argument('-s', '--selection', required=True,
                        help='The specific classification to search for for read selection')
    parser.add_argument('-a', '--all_subtaxa', action='store_true',
                        required=False,
                        help='Use this flag to extract all reads classified as taxa which\n'
                        'fall under the provided classification')
    parser.add_argument('-k', '--kraken_report', required=False,
                        help='The file generated by Kraken2 report flag. Only needed\n'
                             'if using the --all_subtaxa flag')
    parser.add_argument('-t', '--threads', default=1,
                        required=False,
                        help='The number of threads to use')
    parser.add_argument('-o', '--output', required=True,
                        help='The path to the directory to hold the output')

    arguments = parser.parse_args()

ExtractReads(arguments)
