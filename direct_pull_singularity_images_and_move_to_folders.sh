singularity pull flag_cbbasic.image docker://ghcr.io/formbio/flag_cbbasic:latest
mv flag_cbbasic.image containers/cbbasic/
cd containers/combinefilter/
singularity build flag_combinefilter.image Singularity.def
cd ../../
singularity pull flag_entap.image docker://ghcr.io/formbio/flag_entap:latest
mv flag_entap.image containers/entap/
singularity pull flag_exonerate.image docker://ghcr.io/formbio/flag_exonerate:latest
mv flag_exonerate.image containers/exonerate/
singularity pull flag_helixer.image docker://ghcr.io/formbio/flag_helixer:latest
mv flag_helixer.image containers/helixer/
singularity pull flag_helixercpu.image docker://ghcr.io/formbio/flag_helixercpu:latest
mv flag_helixercpu.image containers/helixercpu/
singularity pull flag_liftoff.image docker://ghcr.io/formbio/flag_liftoff:latest
mv flag_liftoff.image containers/liftoff/
singularity pull flag_ncbiclibraries.image docker://ghcr.io/formbio/flag_ncbiclibraries:latest
mv flag_ncbiclibraries.image containers/ncbiclibraries/
singularity pull flag_ncbitools.image docker://ghcr.io/formbio/flag_ncbitools:latest
mv flag_ncbitools.image containers/ncbitools/
singularity pull flag_pasa.image docker://ghcr.io/formbio/flag_pasa:latest
mv flag_pasa.image containers/pasa/
singularity pull flag_tetools.image docker://ghcr.io/formbio/flag_tetools:latest
mv flag_tetools.image containers/tetools/
singularity pull flag_transdecoder.image docker://ghcr.io/formbio/flag_transdecoder:latest
mv flag_transdecoder.image containers/transdecoder/
singularity pull flag_trinity.image docker://ghcr.io/formbio/flag_trinity:latest
mv flag_trinity.image containers/trinity/
