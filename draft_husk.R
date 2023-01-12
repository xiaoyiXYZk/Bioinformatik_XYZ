#subpop draft


load(url("https://ndownloader.figshare.com/files/25950872"))

husk.data = data_humanSkin$data
husk.meta = data_humanSkin$meta
husk <- CreateSeuratObject(counts = husk.data, meta.data = husk.meta)

head(husk@meta.data, 889)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#husk[["subpopulation"]] <- PercentageFeatureSet(husk, pattern = "^MT-")


#seurat_object <- AddMetaData(seurat_object, "group", c("A", "B"))
#husk <- AddMetaData(seurat_object, "group", c("A", "B"))

# Define the string to search
mystring <- "this is a string with foo and bar in it and fooooo and barrrrrr"

# Define the pattern to search for
pattern <- "foo|bar"

# Use the grep function to extract the matching substrings
matches <- grep(pattern, mystring)

# Print the result
print(matches)

# Use the regmatches function to extract the actual substrings
substrings <- regmatches(mystring, matches)

# Print the result
# print(substrings)

