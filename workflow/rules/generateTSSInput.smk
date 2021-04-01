rule getGeneTSS:
        input:
                geneList = lambda wildcard: config["predDir"]+str(preds_config_file.loc[wildcard.pred,"genes"])
        output:
                geneTSS = expand("{outdir}{{pred}}/geneTSS.500bp.bed", outdir=config["outDir"])
        params:
                chrSizes = config["chrSizes"]
        run:

                shell(
                        """
                        cat {input.geneList} | perl -lane 'if  ($F[5] == "+" ) {{print $F[0]."\t".$F[1]."\t".$F[1]."\t".$F[3]."\t".$F[4]."\t".$F[5]}} else {{print $F[0]."\t".$F[2]."\t".$F[2]."\t".$F[3]."\t".$F[4]."\t".$F[5]}}' > {output.geneTSS}.tmp
                        bedtools slop -b 250 -i {output.geneTSS}.tmp -g {params.chrSizes} > {output.geneTSS}
                        """)


