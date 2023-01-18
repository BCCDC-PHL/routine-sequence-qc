def getFastqSubdir(instrument_type) { 
    if (instrument_type == "miseq") {
        def run_dir = new File(params.run_dir.toString())
	if (!run_dir.exists()) {
           System.out.println("Run dir does not exist: ${params.run_dir}")
           System.exit(1)
	}
        def alignment_dirs = []
        run_dir.eachFile {
	    if (it.isDirectory() && it.getName() ==~ /Alignment_\d+/) {
	        alignment_dirs.add(it.name)
            }
	}
        if (alignment_dirs.size() > 0) {
	    def alignment_dir = new File(run_dir.absolutePath + '/' + alignment_dirs.last())
            def timestamp_dirs = []
            alignment_dir.eachFile{
	        if (it.isDirectory() && it.getName() ==~ /\d+_\d+/) {
		    timestamp_dirs.add(it.name)
		}
	    }
            def timestamp_dir = timestamp_dirs.last()
            fastq_subdir = '/' + alignment_dir.name + '/' + timestamp_dir + '/Fastq'
	    println(fastq_subdir)
	    System.exit(0)
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


def makeFastqSearchPath (illumina_suffixes, fastq_exts, instrument_type) {
    def fastq_search_path = []
    def fastq_subdir = getFastqSubdir(instrument_type)
    for (suffix in illumina_suffixes){
        for (ext in fastq_exts){
            fastq_search_path.add(params.run_dir.toString() + fastq_subdir + '/' + suffix.toString() + ext.toString())
        }
    }
    return fastq_search_path
}