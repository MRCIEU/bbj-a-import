library(dplyr)
library(tidyr)
library(jsonlite)
library(parallel)

combine_dat <- function(file1, file2, file3)
{
	cmd <- paste0("cat ", file1, " <( zcat ", file2, " | sed 1d | gzip -c ) > ", file3)
	print(cmd)
	system2("/bin/bash", args = c("-c", shQuote(cmd)))
}

combine_chr <- function(file1, file2)
{
	ft <- paste0(file1, "temp")
	dir.create(ft)
	system(paste0("tar -xvf ", file1, " -C ", ft))

}

config <- read_json('config.json')

a <- read.csv("jbb-a-metadata.csv", stringsAsFactors=FALSE)
a <- subset(a, !is.na(X...newid))
a <- subset(a, !grepl("sex stratified", Study))

a$filename <- ""
a$filename[a$X...newid <= 70] <- paste0("qtl", a$ID[a$X...newid <= 70], ".txt.gz")
a$filename[a$X...newid > 70] <- paste0("cc", a$ID[a$X...newid > 70], ".txt.gz")

a$gwasid <- paste0("jbb-a-", a$X...newid)
a$sample.size <- a$sample.size %>% gsub(",","",.) %>% as.numeric()
a$case <- a$case %>% strsplit(., " ") %>% sapply(., function(x) x[1]) %>% gsub(",","",.) %>% as.numeric
a$ctrl <- a$ctrl %>% strsplit(., " ") %>% sapply(., function(x) x[1]) %>% gsub(",","",.) %>% as.numeric
a$sample.size[is.na(a$sample.size)] <- a$ctrl[is.na(a$sample.size)] + a$case[is.na(a$sample.size)]
a$sex  <-  "Males and females"
a$sex[grepl("female", a$Study)] <- "Females"
a$sex[grepl(", male", a$Study)] <- "Males"
a$category <- "Continuous"
a$category[!is.na(a$case)] <- "Binary"
a$subcategory <- NA
a$unit <- NA
a$group_name <- "public"
a$access <- "public"
a$build <- "HG19/GRCh37"
a$author <- "Ishigaki K"
a$year <- 2019
a$population <- "East Asian"
a$population[grepl("European", a$Study)] <- "European"
a$population[grepl("TransEthnic", a$Study)] <- "Mixed"
a$trait <- a$Study %>% strsplit(., "\\(") %>% sapply(., function(x) x[1]) %>% trimws()
a$trait[grepl("glutamyl trans", a$trait)] <- "Gamma glutamyl transferase"
a$pmid <- a$PMID

tocombine <- a$gwasid[duplicated(a$gwasid)]
for(i in 1:length(tocombine))
{
	message(i)
	x <- subset(a, gwasid == tocombine[i])
	print(x$Study)
	combine_dat(
		file.path(config$datadir, "dl", x$filename[1]), 
		file.path(config$datadir, "dl", x$filename[2]),
		file.path(config$datadir, "dl", paste0(tocombine[i], ".txt.gz"))
	)
	a$filename[a$gwasid == tocombine[i]] <- paste0(tocombine[i], ".txt.gz")
}

a <- subset(a, !duplicated(a$gwasid))

a[["delimiter"]] <- "tab"
a[["header"]] <- TRUE
a[["mr"]] <- 1

a[["chr_col"]] <- 1
a[["pos_col"]] <- 2
a[["oa_col"]] <- 4
a[["ea_col"]] <- 3
a[["snp_col"]] <- 0
a[["eaf_col"]] <- 5
a[["beta_col"]] <- 9
a[["se_col"]] <- 11
a[["pval_col"]] <- 13

index <- a$category == "Continuous"
a[["chr_col"]][index] <- 1
a[["pos_col"]][index] <- 2
a[["oa_col"]][index] <- 3
a[["ea_col"]][index] <- 4
a[["snp_col"]][index] <- 0
a[["eaf_col"]][index] <- 5
a[["beta_col"]][index] <- 7
a[["se_col"]][index] <- 8
a[["pval_col"]][index] <- 9
a$id <- a$gwasid
a$sample_size <- a$sample.size
a$ncase <- a$case
a$ncontrol <- a$ctrl

b <- subset(a, select=c(chr_col, pos_col, oa_col, ea_col, snp_col, eaf_col, beta_col, se_col, pval_col, sex, category, subcategory, unit, population, group_name, access, build, author, year, trait, pmid, id, sample_size, ncase, ncontrol, filename, header, delimiter, mr))

dir.create(file.path(config$datadir, "ready"))
dir.create(file.path(config$datadir, "processed"))
dir.create("job_reports")

for(i in 1:nrow(b))
{
	message(i)
	file.copy(
		file.path(config$datadir, "dl", b$file[i]),
		file.path(config$datadir, "ready", b$file[i])
	)
}

temp <- mclapply(1:nrow(b), function(i)
{
	message(i)
	cmd <- paste0("zcat ", file.path(config$datadir, "ready", b$file[i]), " | wc -l")
	system2("/bin/bash", args = c("-c", shQuote(cmd)), stdout=TRUE) %>% as.numeric() %>% {. - 1}
}, mc.cores=16) %>% unlist()
b$nsnp <- temp

write.csv(b, file="input.csv")
write.csv(b, file=file.path(config$datadir, "ready", "input.csv"))

