include {
    gather_handover_dataset;
} from "./functions.nf"

workflow GATHER_DATA{
    take:
    main:
        log.info "---Running data gathering and transfer to Web ---"
}