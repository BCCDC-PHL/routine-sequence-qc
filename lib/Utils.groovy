import groovy.io.FileType

class Utils {

    public static String getFastqSubdir(instrument_type) { 
	if (instrument_type == "miseq") {
	    def run_dir = new File(params.run_dir.toString())
	    def alignment_dirs = []
	    run_dir.eachFileMatch(FileType.DIRECTORIES, ~/Alignment_\d+/, { it -> alignment_dirs.add( it.name) })
	    if (alignment_dirs.size() > 0) {
		def alignment_dir = new File(run_dir.absolutePath + '/' + alignment_dirs.last())
		def timestamp_dirs = []
		alignment_dir.eachFileMatch(FileType.DIRECTORIES, ~/\d+_\d+/, { it -> timestamp_dirs.add( it.name) })
		def timestamp_dir = timestamp_dirs.last()
		fastq_subdir = '/' + alignment_dir.name + '/' + timestamp_dir + '/Fastq'
	    } else {
		fastq_subdir = "/Data/Intensities/BaseCalls"
	    }
	} else if (instrument_type == "nextseq") {
	    def analysis_dirs = []
	    def run_analysis_dir = new File(params.run_dir.toString() + "/Analysis")
	    run_analysis_dir.traverse(type: groovy.io.FileType.DIRECTORIES, maxDepth: 0) { analysis_dirs.add(it) }
	    def latest_analysis_dir = analysis_dirs.last()
	    def latest_analysis_dir_number = latest_analysis_dir.getName()
	    fastq_subdir = "/Analysis/" + latest_analysis_dir_number + "/Data/fastq"
	} else {
	    System.out.println("Unsupported instrument type: ${instrument_type}")
	    System.exit(1)
	}
	return fastq_subdir
    }

    public static List makeFastqSearchPath ( illumina_suffixes, fastq_exts, instrument_type ) {
	def fastq_search_path = []
	def fastq_subdir = self.getFastqSubdir(instrument_type)
	for (suffix in illumina_suffixes){
            for (ext in fastq_exts){
		fastq_search_path.add(params.run_dir.toString() + fastq_subdir + '/' + suffix.toString() + ext.toString())
            }
	}
	return fastq_search_path
    }
}
