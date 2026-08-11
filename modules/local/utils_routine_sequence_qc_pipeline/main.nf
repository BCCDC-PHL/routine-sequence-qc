import java.nio.file.*

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
    } else if (instrument_type == "i100") {
        def analysis_dirs = []
        def run_analysis_dir = new File(params.run_dir.toString() + "/Analysis")
        run_analysis_dir.traverse(type: groovy.io.FileType.DIRECTORIES, maxDepth: 0) { analysis_dirs.add(it) }
        def latest_analysis_dir = analysis_dirs.last()
        def latest_analysis_dir_number = latest_analysis_dir.getName()
        fastq_subdir = "/Analysis/" + latest_analysis_dir_number + "/Data/BCLConvert/fastq"
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

def findSampleSheet (instrument_type) {
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
	    def timestamp_subdir = '/' + alignment_dir.name + '/' + timestamp_dir
	    def matches = []
            new File(params.run_dir.toString() + timestamp_subdir).traverse(type: groovy.io.FileType.FILES, maxDepth: 0, nameFilter: ~/SampleSheet.*\.csv/) {
                matches << it
            }
            if (matches.size() > 1) {
                log.warn "Multiple SampleSheet*.csv found. using first: ${matches}"
            }
	    samplesheet_path = matches ? matches[0] : null
        } else {
            samplesheet_path = run_dir.absolutePath + '/' + "SampleSheet.csv"
        }
    } else if (instrument_type == "nextseq") {
        def analysis_dirs = []
        def run_analysis_dir = new File(params.run_dir.toString() + "/Analysis")
        run_analysis_dir.traverse(type: groovy.io.FileType.DIRECTORIES, maxDepth: 0) { analysis_dirs.add(it) }
        def latest_analysis_dir = analysis_dirs.last()
        def latest_analysis_dir_number = latest_analysis_dir.getName()
        data_subdir = "/Analysis/" + latest_analysis_dir_number + "/Data"
	def matches = []
        new File(params.run_dir.toString() + data_subdir).traverse(type: groovy.io.FileType.FILES, maxDepth: 0, nameFilter: ~/SampleSheet.*\.csv/) {
            matches << it
        }
        if (matches.size() > 1) {
            log.warn "Multiple SampleSheet*.csv found. using first: ${matches}"
        }
	samplesheet_path = matches ? matches[0] : null
    } else if (instrument_type == "i100") {
        def analysis_dirs = []
        def run_analysis_dir = new File(params.run_dir.toString() + "/Analysis")
        run_analysis_dir.traverse(type: groovy.io.FileType.DIRECTORIES, maxDepth: 0) { analysis_dirs.add(it) }
        def latest_analysis_dir = analysis_dirs.last()
        def latest_analysis_dir_number = latest_analysis_dir.getName()
        inputs_subdir = "/Analysis/" + latest_analysis_dir_number + "/inputs"
	def matches = []
        new File(params.run_dir.toString() + inputs_subdir).traverse(type: groovy.io.FileType.FILES, maxDepth: 0, nameFilter: ~/SampleSheet.*\.csv/) {
            matches << it
        }
        if (matches.size() > 1) {
            log.warn "Multiple SampleSheet*.csv found. using first: ${matches}"
        }
	samplesheet_path = matches ? matches[0] : null
    } else {
        System.out.println("Unsupported instrument type: ${instrument_type}")
        System.exit(1)
    }
    return samplesheet_path
}