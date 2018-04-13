

INDIR=~/WS/nuc/reads/short
OUTDIR=~/WS/nuc/dup_out

INFILES=`echo $INDIR/*`
echo $INFILES

#mkdir -p $OUTDIR

LOG=$OUTDIR/log.txt

THREADS=16
CORRECTK=21


processprobe () {
	probe=$1
	PREFIX=$probe

	sga preprocess -o $OUTDIR/$PREFIX.preprocessed.fa --suffix ":$probe" $INDIR/${probe}_short.fq

	echo -e '\n' `date` 'preprocessed' $probe '\n'


	PREFIX="$PREFIX.preprocessed"
	sga index -a ropebwt --no-reverse -t ${THREADS} $OUTDIR/$PREFIX.fa

	mv $PREFIX.sai $OUTDIR/$PREFIX.sai
	mv $PREFIX.bwt $OUTDIR/$PREFIX.bwt

	echo -e '\n' `date` 'indexed 1' $probe '\n'


	sga correct -t $THREADS -k $CORRECTK --prefix=$OUTDIR/$PREFIX -o $OUTDIR/$PREFIX.ec.fa $OUTDIR/$PREFIX.fa

	echo -e '\n' `date` 'error corrected' $probe '\n'


	PREFIX=$PREFIX.ec

	sga index -a ropebwt -t ${THREADS} $OUTDIR/$PREFIX.fa

	mv $PREFIX.sai $OUTDIR/
	mv $PREFIX.bwt $OUTDIR/
	mv $PREFIX.rsai $OUTDIR/
	mv $PREFIX.rbwt $OUTDIR/	

	echo -e '\n' `date` 'indexed 2' $probe '\n'


	sga filter --no-duplicate-check -t $THREADS --kmer-size $CORRECTK --prefix=$OUTDIR/$PREFIX $OUTDIR/$PREFIX.fa

	echo -e '\n' `date` 'filtered' $probe '\n'
	

	PREFIX=$PREFIX.filter.pass

	mv $PREFIX.sai $OUTDIR/
	mv $PREFIX.bwt $OUTDIR/
	mv $PREFIX.rsai $OUTDIR/
	mv $PREFIX.rbwt $OUTDIR/

}

for probe in '6685_04-06-2015' '6685_16-06-2015'; do
	processprobe $probe &
done
wait

rmdup () {
	probe=$1
	PREFIX=$OUTDIR/$
}

echo 'done'