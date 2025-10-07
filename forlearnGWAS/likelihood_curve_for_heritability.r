# y is the phenotype vector
# K is kinship matrix. Should be approx 1 on main diag and approx 0 elsewhere
# colanmes and rownames of K must be the sample ids
# names(y) must be sample ids
# the function uses the sample ids to match up everything and remove missing values

# retunrs both the mle for the heritability and its stderr and a vector giving the log likelighood computed across the interval 0,1

estimate.mixed.model <- function( y, kinship, make.positive=TRUE ) {
	#data filter begin
	y = y[!is.na(y)]
	if ( length(y) < 1 ) return(NULL)
	use.ids = intersect( names(y), colnames(kinship))
	match.kinship = match( use.ids, colnames(kinship), nomatch=0)
	K = kinship[match.kinship,match.kinship]

	match.phen = match( use.ids, names(y), nomatch=0)
	y = y[match.phen]
	y = scale(y, center=TRUE, scale = FALSE)
	var_y = var(y)
	#data filter end

	K.eigen.trunc = K.eigen = eigen(K,symmetric=TRUE)
	if ( make.positive ) {
		w.eigen = which( K.eigen$values/K.eigen$values[1] > 1.0e-8 )
		neigen = length(w.eigen) # number of non-trivial principal components
		K.eigen.trunc = K.eigen
		K.eigen.trunc$values = K.eigen.trunc$values[w.eigen]
		K.eigen.trunc$vectors = K.eigen.trunc$vectors[,w.eigen]
	}

	z = t(K.eigen.trunc$vectors) %*% y # transformed phenotype
	zz = z*z
	lambda = K.eigen.trunc$values
	cat("estimating heritability\n")

	opt.theta1 = optim( c( 0.1 ),
		fn =function(theta, zz, lambda ) { # log likelihood for mvn
			u = theta[1]*lambda + 1 -theta[1]
			u = var_y * u
			u = ifelse ( abs(u)<1.0e-10, 1.0e-10 , u )
			return(sum(zz/u) + sum( log(u)))
		},
		gr=function( theta, zz, lambda){ # gradient of log L in terms of genetic and environmental variance
			u = theta[1]*lambda + 1 - theta[1]
			u = var_y * u
			u = ifelse ( abs(u)<1.0e-10 ,1.0e-10, u )
			v = (lambda-1)/u
			g1 = sum( (1-zz/u)*v )
			return(c(g1))
		},
		zz, lambda,
		method="L-BFGS-B", lower = c( 1.0e-6), upper = c( 1.0 ), hessian=TRUE )

	mixed.model.data = list( opt.theta1=opt.theta1, vg=opt.theta1$par[1], vg.se = sqrt(1/opt.theta1$hessian[1,1]) )
	
	mixed.model.data$logL.curve = likelihood.curve( zz, lambda, var_y )
	return(mixed.model.data)
}

likelihood.surface <- function( zz, lambda, steps=200 ) {
	step = 1/steps
	logL = matrix( 0, nrow=steps, ncol=steps )
	for( n in 1:steps ) {
		vg = n/steps
		for( m in 1:steps ) {
			ve = m/steps
			u = vg*lambda+ve
			logL[n,m] = 0.5*(sum(zz/u) + sum( log(u)) + length(zz)*log(2*pi))
		}
	}
	return(logL)
}

likelihood.curve <- function( zz, lambda, var_y, steps=200 ) {
	step = 1/steps
	logL = numeric( length=steps )
	for( n in 1:steps ) {
		vg = n/steps
		ve = 1-vg
		u = vg*lambda+ve
		u = var_y*u
		logL[n] = 0.5*(sum(zz/u) + sum( log(u)) + length(zz)*log(2*pi))
	}
	v = (1:steps)/steps
	df = data.frame( vg=v, logL=logL)
	return(df)
}

# R script to read the GRM binary file
ReadGRMBin=function(prefix, size=4){

	BinFileName=paste(prefix,".grm.bin",sep="")
	IDFileName=paste(prefix,".grm.id",sep="")
	id = read.table(IDFileName)
	n=dim(id)[1]

	BinFile=file(BinFileName, "rb");
	grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)

	m <- matrix(nrow = n, ncol = n)
	
	index = 0
	for (index_a in c(1:n)) {
		for (index_b in c(1:n)){
			if( index_a>= index_b ){
				index = index + 1
				m[index_a, index_b]=grm[index]
			}
		}
	}
	for (index_a in c(1:n)) {
		for (index_b in c(1:n)){
			if( index_a<index_b ){
				m[index_a, index_b]=m[index_b, index_a]
			}
		}
	}
	colnames(m) <- id$V1
	rownames(m) <- id$V1
	return(m)
}

# R script to read the plink matrix file
ReadPlinkiMatrix=function(prefix, prefix2){
	k = read.table(prefix)
	name = read.table(prefix2)
	colnames(k)=name$V1
	rownames(k)=name$V1
	return(k)
}



y <- read.table("/media/bs674/it3/PPIproject/LostPPIdecreaseFitness/282_BLUPs_0607_all_traits_HMP32Names.txt", sep="\t",  header=T)
name=y$TaxaName
y <- y$TotalKernelWeight0607summer_blup
names(y) = name

k = ReadPlinkiMatrix("snp.hBN.kinf", "snp.tfam")
estimate <- estimate.mixed.model(y, k)
png("snp.png")
plot(estimate$logL.curve$vg, estimate$logL.curve$logL, main =estimate$vg )
dev.off()

k = ReadPlinkiMatrix("coding_gene.hBN.kinf", "coding_gene.tfam")
estimate <- estimate.mixed.model(y, k)
png("coding_gene.png")
plot(estimate$logL.curve$vg, estimate$logL.curve$logL, main =estimate$vg )
dev.off()

k = ReadPlinkiMatrix("coding_gene_rare.hBN.kinf", "coding_gene_rare.tfam")
estimate <- estimate.mixed.model(y, k)
png("coding_gene_rare.png")
plot(estimate$logL.curve$vg, estimate$logL.curve$logL, main =estimate$vg )
dev.off()


k = ReadPlinkiMatrix("coding_gene_buried.hBN.kinf", "coding_gene_buried.tfam")
estimate <- estimate.mixed.model(y, k)
png("coding_gene_buried.png")
plot(estimate$logL.curve$vg, estimate$logL.curve$logL, main =estimate$vg )
dev.off()

