//
// This file holds several functions specific to the workflow/hic.nf in the nf-core/hic pipeline
//

import nextflow.Nextflow
import nextflow.Channel
import groovy.text.SimpleTemplateEngine
import groovyx.gpars.dataflow.DataflowWriteChannel

class WorkflowHic {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        // digestion parameters
        if (params.digest && params.digestion && !params.digest.containsKey(params.digestion)) {
            Nextflow.error "Unknown digestion protocol. Currently, the available digestion options are ${params.digest.keySet().join(", ")}. Please set manually the '--restriction_site' and '--ligation_site' parameters."
        }

        checkParamIntList(params.bin_size, log)
        checkParamIntList(params.res_dist_decay, log)
        checkParamIntList(params.res_tads, log)
        checkParamIntList(params.res_compartments, log)

        // Check Digestion or DNase Hi-C mode
        //if (!params.dnase && !params.ligation_site) {
        //  Nextflow.error "Ligation motif not found. Please either use the `--digestion` parameters or specify the `--restriction_site` and `--ligation_site`. For DNase Hi-C, please use '--dnase' option"
        //}

    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

        // Check the params 'list of Integer' or Integer (ex bin_size)
    public static void checkParamIntList(def param2check, log) {
        if (param2check !instanceof Integer && !(param2check instanceof String && param2check ==~ /(\d+)(,\d+)*/)){
            println "\n"
            log.error "ERROR: ${param2check} must be integer or list of integer"
            Nextflow.error('Exiting!')
        }
    }

    // In hic.nf: check if the param is int. If true, don't splitCsv and flatten to avoid error
    public static DataflowWriteChannel checkIfInt(def param2check) {
        if (param2check instanceof Integer) {
            return Channel.from( param2check ).toInteger()
        }
        else {
            return Channel.from( param2check ).splitCsv().flatten().toInteger()
        }
    }
}
