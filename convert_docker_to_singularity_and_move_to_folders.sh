singularity build flag_cbbasic.image docker-daemon://gcr.io/bioinfo-devel/flag_cbbasic:latest
mv flag_cbbasic.image containers/cbbasic/
singularity build flag_combinefilter.image docker-daemon://gcr.io/bioinfo-devel/flag_combinefilter:latest
mv flag_combinefilter.image containers/combinefilter/
singularity build flag_entap.image docker-daemon://gcr.io/bioinfo-devel/flag_entap:latest
mv flag_entap.image containers/entap/
singularity build flag_exonerate.image docker-daemon://gcr.io/bioinfo-devel/flag_exonerate:latest
mv flag_exonerate.image containers/exonerate/
singularity build flag_helixer.image docker-daemon://gcr.io/bioinfo-devel/flag_helixer:latest
mv flag_helixer.image containers/helixer/
singularity build flag_helixercpu.image docker-daemon://gcr.io/bioinfo-devel/flag_helixercpu:latest
mv flag_helixercpu.image containers/helixercpu/
singularity build flag_liftoff.image docker-daemon://gcr.io/bioinfo-devel/flag_liftoff:latest
mv flag_liftoff.image containers/liftoff/
singularity build flag_ncbiclibraries.image docker-daemon://gcr.io/bioinfo-devel/flag_ncbiclibraries:latest
mv flag_ncbiclibraries.image containers/ncbiclibraries/
singularity build flag_ncbitools.image docker-daemon://gcr.io/bioinfo-devel/flag_ncbitools:latest
mv flag_ncbitools.image containers/ncbitools/
singularity build flag_pasa.image docker-daemon://gcr.io/bioinfo-devel/flag_pasa:latest
mv flag_pasa.image containers/pasa/
singularity build flag_tetools.image docker-daemon://gcr.io/bioinfo-devel/flag_tetools:latest
mv flag_tetools.image containers/tetools/
singularity build flag_transdecoder.image docker-daemon://gcr.io/bioinfo-devel/flag_transdecoder:latest
mv flag_transdecoder.image containers/transdecoder/
singularity build flag_trinity.image docker-daemon://gcr.io/bioinfo-devel/flag_trinity:latest
mv flag_trinity.image containers/trinity/
