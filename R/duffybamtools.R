
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  File   : duffybamtools.R                                                         #
#  Date   : 12.Mar.2012                                                         #
#  Content: R-Source for package rbamtools                                      #
#  Version: 2.4.2                                                               #
#  Author : W. Kaisers                                                          #
#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #



.onUnload<-function(libpath) { library.dynam.unload("DuffyNGS",libpath) }

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Declaration of generics
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
setGeneric("filename", function(object) standardGeneric("filename"))
setGeneric("isOpen",function(con,rw="") standardGeneric("isOpen"))
setGeneric("bamClose", function(object) standardGeneric("bamClose"))
setGeneric("bamHeader",function(object)standardGeneric("bamHeader"))
setGeneric("getHeader",function(object)standardGeneric("getHeader"))
setGeneric("getRefCount",function(object) standardGeneric("getRefCount"))
setGeneric("getRefData",function(object) standardGeneric("getRefData"))
setGeneric("getRefCoords",function(object,sn)standardGeneric("getRefCoords"))
setGeneric("getNextChunk",function(object,n=1000,alignedOnly=FALSE,primaryOnly=FALSE)standardGeneric("getNextChunk"))

setGeneric("create.index",function(object,idx_filename) standardGeneric("create.index"))
setGeneric("load.index",function(object,filename) standardGeneric("load.index"))
setGeneric("index.initialized",function(object) standardGeneric("index.initialized"))
setGeneric("bamSort",function(object,prefix,byName=FALSE,maxmem=1e+9) standardGeneric("bamSort"))
setGeneric("bamCopy",function(object,writer,refids,verbose=FALSE)standardGeneric("bamCopy"))
setGeneric("bamSave",function(object,...) standardGeneric("bamSave"))
setGeneric("bamWriter",function(x,filename)standardGeneric("bamWriter"))
setGeneric("rewind",function(object)standardGeneric("rewind"))
setGeneric("removeSeqs",function(x,rows)standardGeneric("removeSeqs"))
setGeneric("addSeq",function(object,SN,LN,AS="",M5=0,SP="",UR="")standardGeneric("addSeq"))
setGeneric("head",function(x,...) standardGeneric("head"))
setGeneric("tail",function(x,...) standardGeneric("tail"))
setGeneric("headerLine",function(object) standardGeneric("headerLine"))
setGeneric("refSeqDict",function(object) standardGeneric("refSeqDict"))
setGeneric("headerReadGroup",function(object)standardGeneric("headerReadGroup"))
setGeneric("headerProgram",function(object)standardGeneric("headerProgram"))
setGeneric("headerLine<-",function(object,value)standardGeneric("headerLine<-"))
setGeneric("refSeqDict<-",function(object,value)standardGeneric("refSeqDict<-"))
setGeneric("headerReadGroup<-",function(object,value)standardGeneric("headerReadGroup<-"))
setGeneric("headerProgram<-",function(object,value)standardGeneric("headerProgram<-"))

setGeneric("getPrevAlign",function(object) standardGeneric("getPrevAlign"))
setGeneric("stepNextAlign",function(object)standardGeneric("stepNextAlign"))
setGeneric("stepPrevAlign",function(object)standardGeneric("stepPrevAlign"))
setGeneric("push_back",function(object,value) standardGeneric("push_back"))
setGeneric("pop_back",function(object) standardGeneric("pop_back"))
setGeneric("push_front",function(object,value) standardGeneric("push_front"))
setGeneric("pop_front",function(object) standardGeneric("pop_front"))
setGeneric("writeCurrentAlign",function(object,value) standardGeneric("writeCurrentAlign"))
setGeneric("insertPastCurrent",function(object,value) standardGeneric("insertPastCurrent"))
setGeneric("insertPreCurrent",function(object,value) standardGeneric("insertPreCurrent"))
setGeneric("moveCurrentAlign",function(object,target) standardGeneric("moveCurrentAlign"))
setGeneric("getNextAlign",function(object, alignedOnly=FALSE,primaryOnly=FALSE) standardGeneric("getNextAlign")) 
setGeneric("as.list",function(x,...) standardGeneric("as.list"))
setGeneric("getHeaderText",function(object,delim="\n") standardGeneric("getHeaderText"))
setGeneric("getVal",function(object,member)standardGeneric("getVal"))
setGeneric("setVal",function(object,members,values)standardGeneric("setVal"))
setGeneric("size",function(object) standardGeneric("size"))
setGeneric("nAligns",   function(object)standardGeneric("nAligns"))
setGeneric("readID",function(object) standardGeneric("readID"))
setGeneric("readLocus",function(object) standardGeneric("readLocus"))
setGeneric("refID",function(object) standardGeneric("refID"))
setGeneric("position",function(object) standardGeneric("position"))
setGeneric("nCigar",function(object) standardGeneric("nCigar"))
setGeneric("cigarData",function(object,readUnits=FALSE) standardGeneric("cigarData"))
setGeneric("mismatchData",function(object,readUnits=FALSE) standardGeneric("mismatchData"))
setGeneric("paired", function(object) standardGeneric("paired"))
setGeneric("mateRefID",function(object) standardGeneric("mateRefID"))
setGeneric("matePosition",function(object) standardGeneric("matePosition"))
setGeneric("insertSize",function(object) standardGeneric("insertSize"))
setGeneric("mapQuality",function(object) standardGeneric("mapQuality"))
setGeneric("alignLength",function(object) standardGeneric("alignLength"))
setGeneric("alignSeq",function(object) standardGeneric("alignSeq"))
setGeneric("alignQual",function(object) standardGeneric("alignQual"))
setGeneric("readSeq",function(object) standardGeneric("readSeq"))
setGeneric("readQual",function(object) standardGeneric("readQual"))
setGeneric("readPhredScores",function(object) standardGeneric("readPhredScores"))
setGeneric("nMismatch",function(object) standardGeneric("nMismatch"))
setGeneric("getTag",function(object,tag) standardGeneric("getTag"))
setGeneric("getAllTags",function(object,sep="::") standardGeneric("getAllTags"))
setGeneric("setTag",function(object,tag, value) standardGeneric("setTag"))
setGeneric("reverseStrand", function(object) standardGeneric("reverseStrand"))
setGeneric("firstInPair", function(object) standardGeneric("firstInPair"))
setGeneric("secondInPair", function(object) standardGeneric("secondInPair"))
setGeneric("cigarString", function(object) standardGeneric("cigarString"))
setGeneric("unmapped", function(object) standardGeneric("unmapped"))
setGeneric("secondaryAlign", function(object) standardGeneric("secondaryAlign"))
setGeneric("mateReverseStrand", function(object) standardGeneric("mateReverseStrand"))
setGeneric("mateUnmapped", function(object) standardGeneric("mateUnmapped"))
setGeneric("modifyAlign",function(object,refID=NULL,pos=NULL,mateRefID=NULL, matePos=NULL, 
		readID=NULL, qseq=NULL, qual=NULL) standardGeneric("modifyAlign"))


###################################################################################################
#                                                                                                 #
# bamReader                                                                                       #
#                                                                                                 #
###################################################################################################

setClass("bamReader",representation(filename="character",reader="externalptr",
				index="externalptr",startpos="numeric"),
         validity=function(object){return(ifelse(is.null(object@reader),FALSE,TRUE))})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Opening and closing a BAM-File for reading
#  
setMethod(f="initialize", signature="bamReader",
          definition=function(.Object,filename){
            .Object@filename<-filename
            .Object@reader=.Call("bam_reader_open",path.expand(filename),PACKAGE="DuffyNGS")
            .Object@startpos=.Call("bam_reader_tell",.Object@reader,PACKAGE="DuffyNGS")
            return(.Object)
          }
)

bamReader<-function(filename,indexname,idx=FALSE,verbose=0){
  if(!is.logical(idx))
    stop("[bamReader] idx must be logical!")
  if(length(idx)>1)
    stop("[bamReader] length(idx) must be 1!")
  if(!is.numeric(verbose))
    stop("[bamReader] verbose must be numeric!")
    
  reader<-new("bamReader",filename)
  if((!idx) && missing(indexname))
  {
    if(verbose[1]==1)
      cat("[bamReader] Opened file '",basename(filename),"'.\n",sep="")
    else if(verbose[1]==2)
      cat("[bamReader] Opened file '",filename,"'.\n",sep="")
    return(reader)
  }
  
  # use indexname or set default
  if(missing(indexname))
    idxfile<-paste(filename,"bai",sep=".")
  else
    idxfile<-indexname   
  load.index(reader,idxfile)
 
  if(verbose[1]==1)
    cat("[bamReader] Opened file '",basename(filename),"' and index '",basename(idxfile),"'.\n",sep="")
  else if(verbose[1]==2)
    cat("[bamReader] Opened file '",filename,"' and index '",idxfile,"'.\n",sep="")
  
  return(reader)
}

setMethod("filename", "bamReader", function(object) return(object@filename))
setMethod("isOpen",signature="bamReader",definition=function(con,rw="")
{ return(!(.Call("is_nil_externalptr",con@reader,PACKAGE="DuffyNGS"))) })
setMethod(f="bamClose",signature="bamReader",definition=function(object) {
  if(!.Call("is_nil_externalptr",object@index,PACKAGE="DuffyNGS"))
  {.Call("bam_reader_unload_index",object@index,PACKAGE="DuffyNGS")}
  invisible(.Call("bam_reader_close",object@reader,PACKAGE="DuffyNGS"))
})

setMethod(f="getNextChunk",signature="bamReader",definition=function(object, n=1000, alignedOnly=FALSE,
		primaryOnly=FALSE) {
	rangeObject <- bamRange( NULL, NULL, FALSE);

  	rangeObject@range <- .Call("bam_reader_get_next_chunk",object@reader, as.integer(n), 
					as.logical(alignedOnly), as.logical(primaryOnly), PACKAGE="DuffyNGS");
	return( rangeObject)
})


#  End: Opening and closing a BAM-File for reading
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Header related functions

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  This is one standard Method for creation of bamHeader                        #
#  and is used as a simple way to pass a header to a new                        #
#  instance of bamWriter                                                        #
setMethod(f="getHeader",signature="bamReader",definition=function(object){
  if(!isOpen(object))
    stop("[getHeader.bamReader] reader must be opened! Check with 'isOpen(reader)'!")
  return(new("bamHeader",.Call("bam_reader_get_header",object@reader))) })
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setMethod(f="getHeaderText",signature="bamReader",definition=function(object){
  if(!isOpen(object))
    stop("[getHeaderText.bamReader] reader must be opened! Check with 'isOpen(reader)'!")
  return(new("bamHeaderText",
             .Call("bam_reader_get_header_text",object@reader,PACKAGE="DuffyNGS")))
})

# getRefCount
setMethod(f="getRefCount",signature="bamReader",definition=function(object) {
  if(!isOpen(object))
    stop("[getRefCount.bamReader] reader must be opened! Check with 'isOpen(reader)'!")
  return(.Call("bam_reader_get_ref_count",object@reader,PACKAGE="DuffyNGS"))})

# getRefData
setMethod(f="getRefData",signature="bamReader",definition=function(object) {
  if(!isOpen(object))
    stop("[getRefData.bamReader] reader must be opened! Check with 'isOpen(reader)'!")
  return(.Call("bam_reader_get_ref_data",object@reader,PACKAGE="DuffyNGS"))})

# getRefCoords: Helper function that returns coordinates of entire ref
# for usage with bamRange, gapList or siteList function
setMethod(f="getRefCoords",signature="bamReader",definition=function(object,sn){
  if(!is.character(sn))
    stop("[getRefCoords] sn must be character!")
  if(length(sn)>1)
    stop("[getRefCoords] sn must have length 1!")
  ref<-getRefData(object)
  id<-which(sn==ref$SN)
  if(length(id)==0)
    stop("[getRefCoords] No match for sn in ref-data-SN!")
  coords<-c(ref$ID[id],0,ref$LN[id])
  names(coords)<-c("refid","start","stop")
  return(c(ref$ID[id],0,ref$LN[id]))
})


#  End Header related functions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Index related functions

# create.index
setMethod(f="create.index",signature="bamReader",definition=function(object,idx_filename)
{
  if(missing(idx_filename))
    idx_filename<-paste(object@filename,".bai",sep="")
  invisible(.Call("bam_reader_create_index",path.expand(object@filename),
                 path.expand(idx_filename),PACKAGE="DuffyNGS"))})

setMethod("load.index",signature="bamReader",definition=function(object,filename){
  if(!is.character(filename))
    stop("[bamReader.load.index] Filename must be character!\n")
  if(!file.exists(filename))
    stop("[bamReader.load.index] Index file \"",filename,"\" does not exist!\n")
  
  # Set index Variable in given bamReader object:
  # Read object name, create expression string and evaluate in parent frame
  reader<-deparse(substitute(object))
  extxt<-paste(reader,"@index<-.Call(\"bam_reader_load_index\",\"",
               path.expand(filename),"\",PACKAGE=\"DuffyNGS\")",sep="")
  eval.parent(parse(text=extxt))
  
  # Return true if bamReader@index!=NULL (parent frame)
  extxt<-paste(".Call(\"is_nil_externalptr\",",reader,"@index,PACKAGE=\"DuffyNGS\")",sep="")
  invisible(!eval.parent(parse(text=extxt)))
})
setMethod("index.initialized", signature="bamReader",definition=function(object)
{return(!(.Call("is_nil_externalptr",object@index,PACKAGE="DuffyNGS")))})

setMethod(f="bamSort",signature="bamReader",
          definition=function(object,prefix="sorted",byName=FALSE,maxmem=1e+9)
          {
            maxmem<-floor(maxmem)
            cat("[bamSort] Filename: ",object@filename,"\n")
            cat("[bamSort] Prefix  : ",prefix,"\n")
            cat("[bamSort] Maxmem  : ",maxmem,"\n")
            cat("[bamSort] By Name : ",byName,"\n")
            .Call("bam_reader_sort_file",object@filename,prefix,maxmem,byName,PACKAGE="DuffyNGS")
            cat("[bamSort] Sorting finished.\n")
            return(invisible(paste(prefix,".bam",sep="")))
          })

#  End Index related functions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# getNextAlign
setMethod(f="getNextAlign",signature="bamReader",definition=function(object, alignedOnly=FALSE, primaryOnly=FALSE)
{
  ans<-.Call("bam_reader_get_next_align",object@reader, as.logical(alignedOnly), as.logical(primaryOnly), PACKAGE="DuffyNGS")
  if(is.null(ans))
    return(invisible(NULL))
  else
    return(new("bamAlign",ans))
})



setMethod("rewind","bamReader",function(object) {return(invisible(.Call("bam_reader_seek",object@reader,object@startpos,PACKAGE="DuffyNGS")))})

setMethod("bamSave","bamReader",function(object,writer){
  if(!is(writer,"bamWriter"))
    stop("[bamSave.bamReader] writer must be 'bamWriter'!")
  if(!isOpen(object))
    stop("[bamSave.bamReader] reader is not open! Check 'isOpen'!")
  if(!isOpen(writer))
    stop("[bamSave.bamReader] writer is not open! Check 'isOpen'!")
  
  # Saving old reading position
  oldpos<-.Call("bam_reader_tell",object@reader,PACKAGE="DuffyNGS")
  # Reset reader to start position
  .Call("bam_reader_seek",object@reader,object@startpos,PACKAGE="DuffyNGS")
  
  nAligns<-.Call("bam_reader_save_aligns",object@reader,writer@writer,PACKAGE="DuffyNGS")
  bm<-Sys.localeconv()[7]
  cat("[bamSave.bamReader] Saving ",format(nAligns,big.mark=bm)," to file '",basename(writer@filename),"' finished.\n",sep="")
  .Call("bam_reader_seek",object@reader,oldpos,PACKAGE="DuffyNGS")
  return(invisible(nAligns))
})

setMethod("bamCopy","bamReader",function(object,writer,refids,verbose=FALSE)
{
  if(!is(writer,"bamWriter"))
    stop("[bamCopy.bamReader] writer must be 'bamWriter'!")
  if(!isOpen(object))
    stop("[bamCopy.bamReader] reader is not open! Check 'isOpen'!")
  if(!isOpen(writer))
    stop("[bamCopy.bamReader] writer is not open! Check 'isOpen'!")
  if(!index.initialized(object))
    stop("[bamCopy.bamReader] reader must have initialized index! Check 'index.initialized'!")
  
  # Check refids argument: When missing copy all ref's
  ref<-getRefData(object)
  if(missing(refids))
  {
    refids<-ref$ID
    n<-length(refids)
    mtc<-1:n
  }
  else
  {
    mtc<-match(refids,ref$ID)
    if(any(is.na(mtc)))
      stop("[bamCopy.bamReader] refids must be subset of Reference-ID's! Check 'getRefData'!")
    n<-length(refids)    
  }

  # Copy aligns with bamRanges as intermediate buffer
  bm<-Sys.localeconv()[7]
  nAligns<-0
  for(i in 1:n)
  {
    range<-bamRange(object,c(ref$ID[mtc[i]],0,ref$LN[mtc[i]]),complex=FALSE)
    nAligns<-nAligns+size(range)
    if(verbose)
      cat("[bamCopy.bamReader] i: ",i,"\tCopying ",format(size(range),big.mark=bm,width=10)," aligns for Reference '",ref$SN[mtc[i]],"'.\n",sep="")
    bamSave(writer,range)
    rm(range)
    gc()    
  }
  cat("[bamCopy.bamReader] Copying ",format(nAligns,big.mark=bm,width=10)," aligns finished.\n",sep="")
})

###################################################################################################
#                                                                                                 #
# bamHeader                                                                                       #
# Description: See SAM File Format Specification (v1.4-r985) September 7,2011, Section 1.3        #
#                                                                                                 #
###################################################################################################

setClass("bamHeader",representation(header="externalptr"),
         validity=function(object){return(ifelse(is.null(object@header),FALSE,TRUE))})

setMethod("initialize","bamHeader",function(.Object,extptr){
  if(!is(extptr,"externalptr"))
    stop("[initialize.bamHeader] extptr must be externalptr!")
  .Object@header<-extptr
  return(.Object)
})

setMethod(f="getHeaderText",signature="bamHeader",definition=function(object) {
  return(new("bamHeaderText",.Call("bam_header_get_header_text",object@header,PACKAGE="DuffyNGS"))) })

setMethod("as.character","bamHeader",function(x,...){
  .Call("bam_header_get_header_text",x@header,PACKAGE="DuffyNGS")
})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  This is the main function for creating an instance of bamWriter              #

setMethod("bamWriter","bamHeader",function(x,filename){
  if(!is.character(filename))
    stop("[bamWriter.bamHeader] filename must be character!")
  return(new("bamWriter",x,filename))
})
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# headerLine: Represents two entries: Format version (VN) and sorting order (SO)
# Valid format for VN : /^[0-9]+\.[0-9]+$/.
# Valid entries for SO: unknown (default), unsorted, queryname, coordinate.

setClass("headerLine",representation(VN="character",SO="character"),
         validity=function(object)
         {
           if(length(VN)==1 & length(SO)==1)
             return(TRUE)
           else
             return(FALSE)
         })

setMethod(f="initialize",signature="headerLine",definition=function(.Object,hl="",delim="\t"){
  # Parses header line from header section
  if(!is.character(hl))
    stop("[headerLine.initialize] Argument must be string.\n")
  
  # Default object content (hl="" or character(0))
  if((length(hl)==1 && nchar(hl)==0) || length(hl)==0)
  {
    .Object@VN="1.4"
    .Object@SO="unknown"
    return(.Object)
  }
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #  Split input string into tags
  tags<-unlist(strsplit(hl,delim))
  #  Three tags!
  if(length(tags)!=3)
    stop("[headerLine.initialize] hl must contain three tags separated by '",delim,"'!\n")
  #  First  tag: '@HD'
  if(tags[1]!="@HD")
    stop("[headerLine.initialize] First tag of string must be @HD!\n")
  #  Second tag: 'VN'
  #  TODO: Check Accepted format: /^[0-9]+\.[0-9]+$/.  
  if(substr(tags[2],1,2)!="VN")
    stop("[headerLine.initialize] Second tag of string must be VN!\n")
  .Object@VN=substring(tags[2],4)
  # Third   tag: 'SO'
  if(substr(tags[3],1,2)!="SO")
    stop("[headerLine.initialize] Third tag of string must be SO!\n")
  
  str<-substring(tags[3],4)
  if(str=="coordinate")
    .Object@SO<-"coordinate"
  else if(str=="unknown")
    .Object@SO<-"unknown"
  else if(str=="unsorted")
    .Object@SO<-"unsorted"
  else if(str=="queryname")
    .Object@SO<-"queryname"
  
  return(.Object)
})

setMethod("getHeaderText","headerLine",function(object,delim="\t")
{return(paste("@HD\tVN:",object@VN,"\tSO:",object@SO,sep=""))})

setMethod("getVal",signature="headerLine",definition=function(object,member){
  if(!is.character(member))
    stop("[getVal.headerLine] Member must be character!\n")
  if(member=="VN")
    return(object@VN)
  if(member=="SO")
    return(object@SO)
  stop("[getVal.headerLine] Member '",member,"' must be 'VN' or 'SO'!\n")
})

setMethod("setVal",signature="headerLine",definition=function(object,members,values){
  if(!is.character(members) || !is.character(values))
    stop("[setVal.headerLine] Members and values must be character!\n")
  if(length(members)!=length(values))
    stop("[setVal.headerLine] Members and values must have same length!\n")
  tagLabs<-c("VN","SO")
  mtc<-match(members,tagLabs)
  if(any(is.na(mtc)))
    stop("[setVal.headerLine] Member names must be valid Header line entries!\n")
  n<-length(members)
  if(n>2)
    stop("[setVal.headerLine] Only two members can be set!\n")
  obj<-deparse(substitute(object))
  for(i in 1:n)
  {
    txt<-paste(obj,"@",members[i],"<-'",values[i],"'",sep="")
    eval.parent(parse(text=txt))
  }
  return(invisible())
})

setMethod("as.list",signature="headerLine",definition=function(x,...)
{return(list(VN=x@VN,SO=x@SO))})

setMethod("show","headerLine",function(object)
{
  cat("An object of class \"",class(object),"\"\n",sep="")
  cat("VN: ",object@VN,"\nSO: ",object@SO,"\n",sep="")
})

#  End headerLine
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  refSeqDict: Reference Sequence Dictionary                                    #
#  Represents a variable number of Ref Seqs                                     #
#  Valid Members (Entries for each sequence, stored in a data.frame):           #
#  SN Reference sequence name                                                   #
#  LN Reference sequence length                                                 #
#  AS Genome assembly identifier                                                #
#  M5 MD5 checksum of the sequence                                              #
#  SP Species                                                                   #
#  UR URI of the sequence                                                       #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("refSeqDict",representation(SN="character",LN="numeric",AS="character",M5="numeric",SP="character",UR="character"))
setMethod(f="initialize",signature="refSeqDict",definition=function(.Object,hsq="",delim="\t")
{
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  # Parses Reference sequence dictionary of header-text
  # hsq= Vector of characters, each representing one Ref-Sequence
  # length(hsq) = number of Ref-Sequences
  # Each Ref-string contains 'internally' [tab] delimited seqments:
  #               "SN:ab\tLN:12\tAS:ab\tM5:12\tSP:ab\tUR:ab"
  # It's allowed to skip segments
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  
  if(!is.character(hsq))
    stop("[refSeqDict.initialize] hsq must be character!")
  
  n<-length(hsq)
  # Return empty object when no input string is given
  if((n==1 && nchar(hsq)==0) || n==0)
    return(.Object)
  
  .Object@SN<-character(n)
  .Object@LN<-numeric(n)
  .Object@AS<-character(n)
  .Object@M5<-numeric(n)
  .Object@SP<-character(n)
  .Object@UR<-character(n)
  labels<-c("SN","LN","AS","M5","SP","UR")
  for(i in 1:n)
  {
    # Containes separated tags for one sequence
    seq<-unlist(strsplit(hsq[i],delim))
    if(seq[1]!="@SQ")
      stop("[initialize.refSeqDict] First segment in Ref-sequence tag must be '@SQ'!")
    seq<-seq[-1]
    
    # Contains column number in dict@df for each tag
    cols<-match(substr(seq,1,2),labels)
    m<-length(cols)
    for(j in 1:m)
    {
      txt<-substr(seq[j],4,nchar(seq[j]))
      # Empty entries are skipped (to avoid errors)
      if(nchar(txt)>0)
      {
        if(cols[j]==1)
          .Object@SN[i]<-txt
        else if(cols[j]==2)
        {
          # Try to convert into numeric value
          numb<-suppressWarnings(as.numeric(txt))  
          if(is.na(numb))
          {warning("[refSeqDict.initialize] No numeric value for LN: '",txt,"'!\n",sep="")}
          else
            .Object@LN[i]<-numb
        }
        else if(cols[j]==3)
          .Object@AS[i]<-txt
        else if(cols[j]==4)
        {
          # Try to convert into numeric value
          numb<-suppressWarnings(as.numeric(txt))  
          if(is.na(numb))
          {warning("[refSeqDict.initialize] No numeric value for LN: '",txt,"'!\n",sep="")}
          else
            .Object@M5<-numb
        }
        else if(cols[j]==5)
          .Object@SP[i]<-txt
        else if(cols[j]==6)
          .Object@UR[i]<-txt
      }
    }
  }
  return(.Object)
})

setMethod(f= "[",signature="refSeqDict",definition=function(x,i){
  rsd<-new("refSeqDict")
  rsd@SN<-x@SN[i]
  rsd@LN<-x@LN[i]
  rsd@AS<-x@AS[i]
  rsd@M5<-x@M5[i]
  rsd@SP<-x@SP[i]
  rsd@UR<-x@UR[i]
  return(rsd)
})

setMethod(f="dim",signature="refSeqDict",definition=function(x){return(c(length(x@SN),6))})

setMethod("removeSeqs",signature="refSeqDict",definition=function(x,rows){
  # Removes given rows (=Sequences) from Dictionary so they are excluded from header
  n<-length(x@SN)
  if(!is.numeric(rows))  
    stop("[removeSeqs.refSeqDict] Sequence indices must be numeric!")
  rows<-as.integer(rows)
  if(any(rows)<1)
    stop("[removeSeqs.refSeqDict] Sequence indices must be positive!")
  if(any(rows)>n)
    stop("[removeSeqs.refSeqDict] Sequence indices must be <",n,"!")
  
  # Execute per eval in parent.frame
  if(length(rows)>1)
    rmv<-paste("c(",paste(rows,collapse=","),")",sep="")
  else
    rmv<-rows
  obj<-deparse(substitute(x))
  dictcol<-paste(obj,"@SN",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@LN",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@AS",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@M5",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@SP",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  
  dictcol<-paste(obj,"@UR",sep="")
  eval.parent(parse(text=paste(dictcol,"<-",dictcol,"[-",rmv,"]",sep="")))
  return(invisible())
})

setMethod("addSeq",signature="refSeqDict",definition=function(object,SN,LN,AS="",M5=0,SP="",UR=""){
  index<-length(object@SN)+1
  obj<-deparse(substitute(object))
  colidx<-paste("[",index,"]",sep="")
  
  # Appends new Sequence (row) at the end
  dictcol<-paste(obj,"@SN",colidx,sep="")
  eval.parent(parse(text=paste(dictcol,"<-'",SN,"'",sep="")))
  
  dictcol<-paste(obj,"@LN",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-",LN,sep="")))
  
  dictcol<-paste(obj,"@AS",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-'",AS,"'",sep="")))
  
  dictcol<-paste(obj,"@M5",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-",M5,sep="")))
  
  dictcol<-paste(obj,"@SP",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-'",SP,"'",sep="")))
  
  dictcol<-paste(obj,"@UR",colidx,sep="") 
  eval.parent(parse(text=paste(dictcol,"<-'",UR,"'",sep="")))
  
  return(invisible())
})

setMethod("getHeaderText",signature="refSeqDict",definition=function(object,delim="\t"){ 
  # Returns Ref Data String (can be used for creating new BAM file via bamWriter)
  labels<-c("SN","LN","AS","M5","SP","UR")
  n<-length(object@SN)
  if(n==0)
    return(character(0))
  
  seqs<-character(n)
  
  for(i in 1:n)
  {
    ans<-"@SQ"    
    if(nchar(object@SN[i])>0)
      ans<-paste(ans,delim,"SN:",object@SN[i],sep="")
    if(object@LN[i]>0)
      ans<-paste(ans,delim,"LN:",object@LN[i],sep="")
    if(nchar(object@AS[i])>0)
      ans<-paste(ans,delim,"AS:",object@AS[i],sep="")
    if(object@M5[i]>0)
      ans<-paste(ans,delim,"M5:",object@M5[i],sep="")
    if(nchar(object@SP[i])>0)
      ans<-paste(ans,delim,"SP:",object@SP[i],sep="")
    if(nchar(object@UR[i])>0)
      ans<-paste(ans,delim,"UR:",object@UR[i],sep="")
    seqs[i]<-ans
  }
  return(paste(seqs,collapse="\n"))
})

# Return first or last part of refSeqDict data.frame
# S3 Generic is supplied via importFrom in NAMESPACE
setMethod("head","refSeqDict",function(x,n=6L,...) {
  stopifnot(length(n) == 1L)
  if (n < 0L)
    stop("[head.refSeqDict] n<0!")
  m<-length(x@SN)
  if(m==0)
    cat("[head.refSeqDict] Empty object.\n")
  
  n<-min(n,m)
  if(n == 0L)
    return(as.data.frame(new("refSeqDict")))
  else
    return(as.data.frame(x)[1:n,])
})
# S3 Generic is supplied via importFrom in NAMESPACE
setMethod("tail","refSeqDict",definition=function(x,n=6L,...) {
  stopifnot(length(n) == 1L)
  if (n < 0L)
    stop("[tail.refSeqDict] n<0!")
  m<-length(x@SN)
  if(m==0)
    cat("[tail.refSeqDict] Empty object.\n") 
  n<-min(n,m)
  if(n == 0L)
    return(as.data.frame(new("refSeqDict")))
  else
  {
    n<-m-n+1
    return(x@df[n:m,])
  }
})

setMethod("show","refSeqDict",function(object){
  cat("An object of class \"",class(object),"\"\n",sep="")
  print(head(object))
})

#  End refSeqDict
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# headerReadGroup
# ReadGroup
# ID Read Group identifier
# CN Name of sequencing center
# DS Description
# FO Flow order
# KS Nucleotides corresponding to key sequence of each read
# LB Library
# PG Programs used for processing the Read Group
# PI Predicted median insert size
# PL Sequencing Platform:
#    CAPILLARY,LS454,ILLUMINA,SOLID,HELICOS,IONTORRENT or PACBIO
# SM Sample name.


setClass("headerReadGroup",representation(l="list"),validity=function(object) {return(TRUE)})

setMethod(f="initialize",signature="headerReadGroup", definition=function(.Object,hrg="",delim="\t"){
  # Parses Read-Group part of Header data. See Sam Format Specificatioin 1.3 (Header Section)
  .Object@l<-list()
  if(!is.character(hrg))
    stop("[headerReadGroup.initialize] Argument must be string.\n")
  # hrg="" or character(0)
  if((length(hrg)==1 && nchar(hrg)==0) || length(hrg)==0)
    return(.Object)
  # Split string into fields
  tags<-unlist(strsplit(hrg,delim))
  if(tags[1]!="@RG")
    stop("[headerReadGroup.initialize] First item of string must be @RG!\n")
  tags<-tags[-1]
  tagLabs<-c("ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM")
  n<-length(tags)
  for(i in 1:n)
  {
    f<-substr(tags[i],1,2)
    mtc<-match(f,tagLabs)
    if(is.na(mtc))
      stop("[headerReadGroup.initialize] Field identifier '",f,"' not in List!\n")
    .Object@l[[f]]<-substr(tags[i],4,nchar(tags[i]))
  }
  return(.Object)
})

setMethod("getHeaderText",signature="headerReadGroup",definition=function(object,delim="\t") {
  n<-length(object@l)
  if(n==0)
    return(character(0))
  rfstr<-character(n)
  for(i in 1:n)
    rfstr[i]<-paste(names(object@l)[i],object@l[[i]],sep=":")
  return(paste("@RG",paste(rfstr,collapse=delim),sep=delim))
})

setMethod("getVal",signature="headerReadGroup",definition=function(object,member) {
  if(!is.character(member))
    stop("[getVal.headerReadGroup] Member must be character!\n")
  tagLabs<-c("ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM")
  mtc<-match(member[1],tagLabs)
  if(is.na(mtc))
    stop("[getVal.headerReadGroup] Invalid member name!\n")
  return(object@l[[member]])
})

setMethod("setVal",signature="headerReadGroup",definition=function(object,members,values){
  if(!is.character(members) || !is.character(values))
    stop("[setVal.headerReadGroup] Member name and value must be character!\n")
  if(length(members)!=length(values))
    stop("[setVal.headerReadGroup] members and values must have same length!\n")
  tagLabs<-c("ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM")
  mtc<-match(members,tagLabs)
  if(any(is.na(mtc)))
    stop("[setVal.headerReadGroup] Members must be valid Read Group Entries (See SAM Format Specification 1.3!\n")
  n<-length(members)
  obj<-deparse(substitute(object))
  for(i in 1:n)
  {
    txt<-paste(obj,"@l$",members[i],"<-'",values[i],"'",sep="")
    eval.parent(parse(text=txt))
  }
  return(invisible())
})

setMethod("as.list",signature="headerReadGroup",definition=function(x,...){return(x@l)})

#  End headerReadGroup
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# headerProgram
setClass("headerProgram",representation(l="list"),validity=function(object){return(TRUE)})

setMethod(f="initialize",signature="headerProgram",
          definition=function(.Object,hp="",delim="\t")
          {
            # Parses Program part of Header data. See Sam Format Specificatioin 1.3 (Header Section)
            .Object@l<-list()
            if(!is.character(hp))
              stop("[headerProgram.initialize] Argument must be string.\n")
            # hp="" or character(0)
            if((length(hp)==1 && nchar(hp)==0)||length(hp)==0)
              return(.Object)
            # Split string into fields
            tags<-unlist(strsplit(hp,delim))
            if(tags[1]!="@PG")
              stop("[headerProgram.initialize] First item of string must be @PG!\n")
            tags<-tags[-1]
            tagLabs<-c("ID","PN","CL","PP","VN")
            n<-length(tags)
            for(i in 1:n)
            {
              f<-substr(tags[i],1,2)
              mtc<-match(f,tagLabs)
              if(is.na(mtc))
                stop("[heaProgram.initialize] Field identifier '",f,"' not in List!\n")
              .Object@l[[f]]<-substr(tags[i],4,nchar(tags[i]))
            }
            return(.Object)
          })

setMethod("getHeaderText",signature="headerProgram",definition=function(object,delim="\t") {
  n<-length(object@l)
  if(n==0)
    return(character(0))
  
  rfstr<-character(n)
  for(i in 1:n)
    rfstr[i]<-paste(names(object@l)[i],object@l[[i]],sep=":")
  return(paste("@PG",paste(rfstr,collapse=delim),sep=delim))
})

setMethod("getVal",signature="headerProgram",definition=function(object,member) {
  if(!is.character(member))
    stop("[getVal.headerProgram] Member must be character!\n")
  tagLabs<-c("ID","PN","CL","PP","VN")
  mtc<-match(member[1],tagLabs)
  if(is.na(mtc))
    stop("[getVal.headerProgram] Invalid member name!\n")
  return(object@l[[member]])
})

setMethod("setVal",signature="headerProgram",definition=function(object,members,values){
  if(!is.character(members) || !is.character(values))
    stop("[setVal.headerProgram] Member name and value must be character!\n")
  if(length(members)!=length(values))
    stop("[setVal.headerProgram] members and values must have same length!\n")
  tagLabs<-c("ID","PN","CL","PP","VN")
  mtc<-match(members,tagLabs)
  if(any(is.na(mtc)))
    stop("[setVal.headerProgram] Members must be valid Program Entries (See SAM Format Specification 1.3!\n")
  n<-length(members)
  obj<-deparse(substitute(object))
  for(i in 1:n)
  {
    txt<-paste(obj,"@l$",members[i],"<-'",values[i],"'",sep="")
    eval.parent(parse(text=txt))
  }
  return(invisible())
})

setMethod("as.list",signature="headerProgram",definition=function(x,...){return(x@l)})

setMethod("show","headerProgram",function(object)
{
  cat("An object of class \"",class(object),"\"\n",sep="")
  for(i in 1:length(object@l))
  {
    cat(names(object@l)[i],":",object@l[[i]],"\n")
  }
  return(invisible())
})


#  End headerProgram
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  bamHeaderText: Represents and manages textual version of bamHeader           #
#  See SAM Format Specification (v1.4-r985)                                     #
#                                                                               #
#  Contains header Segments :                                                   #
#   head  = headerLine        : @HD Header Line                                 #
#   dict  = refSeqDict        : @SQ Reference Sequence dictionary               #
#   group = headerReadGroup   : @RG Read Group                                  #
#   prog  = headerProgram     : @PG Program                                     #
#                                                                               #
#   TODO:
#   com   = headerComment     : @CO One-line text comment                       #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Class definition and creational routines for bamHeaderText

setClass("bamHeaderText",representation(head="headerLine",dict="refSeqDict",
                                        group="headerReadGroup",prog="headerProgram",com="character"))

setMethod(f="initialize",signature="bamHeaderText", definition=function(.Object,bh="",delim="\n")
{
  # Parses Header data (as reported by getHeaderText)
  # See Sam Format Specification 1.3 (Header Section)
  if(!is.character(bh))
    stop("[bamHeaderText.initialize] Argument must be string.\n")
  
  # Create empty header Set (so it's legal to call getHeaderText()')
  if(length(bh)==1 && nchar(bh)==0)
  {
    .Object@head<-new("headerLine")
    .Object@dict<-new("refSeqDict")
    .Object@group<-new("headerReadGroup")
    .Object@prog<-new("headerProgram")
    return(.Object)
  }
  
  # Split input string: Each fragment contains data for one header segment
  bht<-unlist(strsplit(bh,split=delim))
  
  # Read Header Line
  bhl<-bht[grep("@HD",bht)]
  .Object@head<-new("headerLine",bhl)
  
  # Read Sequence Directory
  bsd<-bht[grep("@SQ",bht)]
  .Object@dict<-new("refSeqDict",bsd)
  
  # Read Group
  brg<-bht[grep("@RG",bht)]
  .Object@group<-new("headerReadGroup",brg)
  
  # Read Program Data
  bpd<-bht[grep("@PG",bht)]
  .Object@prog<-new("headerProgram",bpd)
  
  # Read Text comment
  btc<-bht[grep("@CO",bht)]
  com<-substring(btc,3)
  return(.Object)
})

bamHeaderText<-function(head=NULL,dict=NULL,group=NULL,prog=NULL,com=NULL)
{
  bh<-new("bamHeaderText")
  if(!is.null(head))
  {
    if(is(head,"headerLine"))
      bh@head<-head
    else
      stop("[bamHeaderText] head must be 'headerLine'!")
  }
  if(!is.null(dict))
  {
    if(is(dict,"refSeqDict"))
      bh@dict<-dict
    else
      stop("[bamHeaderText] dict must be 'refSeqDict'")
  }
  if(!is.null(group))
  {
    if(is(group,"headerReadGroup"))
      bh@group<-group
    else
      stop("[bamHeaderText] group must be 'headerReadGroup'!")
  }
  if(!is.null(prog))
  {
    if(is(prog,"headerProgram"))
      bh@prog<-prog
    else
      stop("[bamHeaderText] prog must be 'headerProgram'!")
  }
  if(!is.null(com))
  {
    if(is.character(com))
      bh@com<-com
    else
      stop("[bamHeaderText] com must be 'character'!")
  }
  return(invisible(bh))
}

#  End: Class definition and creational routines for bamHeaderText
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  Public accessors for member objects for bamHeaderText
setMethod(f="headerLine",signature="bamHeaderText",definition=function(object) {return(object@head)})
setMethod(f="refSeqDict",signature="bamHeaderText",definition=function(object) {return(object@dict)})
setMethod(f="headerReadGroup",signature="bamHeaderText",definition=function(object){return(object@group)})
setMethod(f="headerProgram",signature="bamHeaderText",definition=function(object){return(object@prog)})

setReplaceMethod("headerLine","bamHeaderText",function(object,value)
{
  if(!is(value,"headerLine"))
    stop("[headerLine<-.bamHeaderText] value must be 'headerLine'!")
  object@head<-value
  return(object)
})

setReplaceMethod("refSeqDict","bamHeaderText",function(object,value)
{
  if(!is(value,"refSeqDict"))
    stop("[refSeqDict<-.bamHeaderText] value must be 'refSeqDict'!")
  object@dict<-value
  return(object)
})

setReplaceMethod("headerReadGroup","bamHeaderText",function(object,value)
{
  if(!is(value,"headerReadGroup"))
    stop("[headerReadGroup<-.bamHeaderText] value must be 'headerReadGroup'!")
  object@group<-value
  return(object)
})

setReplaceMethod("headerProgram","bamHeaderText",function(object,value)
{
  if(!is(value,"headerProgram"))
    stop("[headerProgram<-.bamHeaderText] value must be 'headerProgram'!")
  object@prog<-value
  return(object)
})

#  End: Public accessors for member objects for bamHeaderText
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #



setMethod("getHeaderText",signature="bamHeaderText",definition=function(object,delim="\n") {
  hd<-getHeaderText(object@head)
  if(length(hd)==0)
    return(character(0))
  hd<-paste(hd,delim,sep="")
  
  dt<-getHeaderText(object@dict)
  if(length(dt)==0)
    return(character(0))
  dt<-paste(dt,delim,sep="")
  
  gp<-getHeaderText(object@group)
  if(length(gp)>0)
    gp<-paste(gp,delim,sep="")
  
  pg<-getHeaderText(object@prog)
  if(length(pg)>0)
    pg<-paste(pg,delim,sep="")
  
  if(length(object@com)>0)
    cm<-paste(paste("@CO",object@com,sep="\t"),collapse=delim)
  else
    cm<-character(0)
  return(paste(hd,dt,gp,pg,cm,sep=""))
})

setMethod("bamHeader","bamHeaderText",
          function(object){return(new("bamHeader",.Call("init_bam_header",getHeaderText(object))))})


###################################################################################################
#                                                                                                 #
# bamWriter class                                                                                 #
# Encapsulates an write-opened Connection to a BAM-file.                                          #
#                                                                                                 #
###################################################################################################

setClass("bamWriter",representation(filename="character",writer="externalptr"),
         validity=function(object) {return(ifelse(is.null(object@writer),FALSE,TRUE))})

setMethod(f="initialize", signature="bamWriter",
          definition=function(.Object,header,filename){
            if(!is(header,"bamHeader"))
              stop("[initialize.bamWriter] header must be bamHeader!\n")
            if(!is.character(filename))
              stop("[initialize.bamWriter] filename must be character!\n")
            .Object@filename<-filename
            .Object@writer<-.Call("bam_writer_open",header@header,filename,PACKAGE="DuffyNGS")
            return(.Object)
          })


setMethod("filename", "bamWriter", function(object) return(object@filename))
setMethod("isOpen",signature="bamWriter",definition=function(con,rw="")
{return(!(.Call("is_nil_externalptr",con@writer,PACKAGE="DuffyNGS")))})

setMethod(f="bamClose",signature="bamWriter",definition=function(object)
{ invisible(.Call("bam_writer_close",object@writer,PACKAGE="DuffyNGS"))})

setMethod(f="bamSave",signature="bamWriter",definition=function(object,value) 
{
  if(is(value,"bamAlign"))
  {
    return(invisible(.Call("bam_writer_save_align",object@writer,value@align,PACKAGE="DuffyNGS")))
  }
  if(is(value,"bamRange"))
  {
    return(invisible(.Call("bam_range_write",object@writer,value@range,PACKAGE="DuffyNGS")))
  }
  else
    stop("bamSave: Saved object must be of type bamAlign or bamRange!\n")
})


###################################################################################################
#                                                                                                 #
# bamRange                                                                                        #
#                                                                                                 #
###################################################################################################

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Encapsulates a bunch of Alignment datasets that typically have been read from a defined         #
# reference region in a BAM-file.                                                                 #
# Technically, the alignments are stored in a (C-implemented) double linked list.                 #
# bamRange objects can be created by a reading procedure on an indexed BAM-file. The alignments   #
# can be iterated, readed, written, deleted and added. bamRange objects can be written to a       #
# BAM-file via an Instance of bamWriter.                                                          #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

bamRange<-function(reader=NULL,coords=NULL,complex=FALSE) {
  if(!is.null(reader))
  {
    if(!index.initialized(reader))
      stop("[bamRange] reader must have initialized index! Use 'load.index'!")
  }
  return(new("bamRange",reader,coords,complex))
}

setClass("bamRange",representation(range="externalptr"),
         validity=function(object) { return(ifelse(is.null(object@range),FALSE,TRUE)) })

setMethod(f="initialize",signature="bamRange",
          definition=function(.Object,reader=NULL,coords=NULL,complex=FALSE)
          {        
            # +++++++++++++++++++++++++++++++++++++++++++            
            #  Create empty range
            if(is.null(reader))
            {
              .Object@range<-.Call("bam_range_init")
              return(.Object)
            }
            
            # +++++++++++++++++++++++++++++++++++++++++++
            #  Create range from bam-file
            if(!is(reader,"bamReader"))
              stop("[bamRange.initialize] reader must be an instance of bamReader!")
            if(length(coords)!=3)
              stop("[bamRange.initialize] coords must be numeric with length=3 (ref,start,stop)!")  
            if(is.null(reader@index))
              stop("[bamRange.initialize] bamReader must have initialized index!")
            if(!is(complex,"logical"))
              stop("[bamRange.initialize] complex must be logical!")
            if(length(complex)>1)
              stop("[bamRange.initialize] complex must have length 1!")
            if(!index.initialized(reader))
              stop("[bamRange.initialize] reader must have initialized index! Use 'load.index'!")
            .Object@range<-.Call("bam_range_fetch",reader@reader,reader@index,trunc(coords),complex,PACKAGE="DuffyNGS")
            return(.Object)
          })

# REMOVED:
# setMethod("as.data.frame",signature="bamRange",definition=function(x,row.names=NULL,optional=FALSE,...) {
#   return(.Call("bam_range_get_align_df",x@range,PACKAGE="DuffyNGS"))
# })



setMethod("size",signature="bamRange",definition=function(object)
{.Call("bam_range_get_size",object@range,PACKAGE="DuffyNGS")})

setMethod("show","bamRange",function(object){
  cat("An object of class '",class(object),"'. size: ",size(object),"\n",sep="")
  return(invisible())
})

setMethod("getNextAlign",signature="bamRange",definition=function(object, alignedOnly=FALSE, primaryOnly=FALSE)
{
  ans<-.Call("bam_range_get_next_align",object@range, as.logical(alignedOnly), as.logical(primaryOnly), PACKAGE="DuffyNGS")
  # Must be checked because align list returns NULL when end is reached
  if(is.null(ans))
    return(ans)
  else
    return(new("bamAlign",ans))
})


setMethod("getPrevAlign",signature="bamRange",definition=function(object)
{ return(new("bamAlign",.Call("bam_range_get_prev_align",object@range,PACKAGE="DuffyNGS")))})

setMethod("stepNextAlign",signature("bamRange"),definition=function(object)
{
  .Call("bam_range_step_next_align",object@range)
  return(invisible())
})

setMethod("stepPrevAlign",signature("bamRange"),definition=function(object)
{
  .Call("bam_range_step_prev_align",object@range)
  return(invisible())
})

# Resets current align to NULL position (i.e. before first element)
# The next call to getNextAlign then returns the first element of list
#setGeneric("windBack", function(object) standardGeneric("windBack"))
setMethod("rewind",signature="bamRange",definition=function(object)
{invisible(.Call("bam_range_wind_back",object@range,PACKAGE="DuffyNGS"))})

setMethod("push_back",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("push_back.bamRange: pushed object must be of class \"bamAlign\"\n")
  .Call("bam_range_push_back",object@range,value@align,PACKAGE="DuffyNGS")
})

setMethod("pop_back",signature="bamRange",definition=function(object)
{.Call("bam_range_pop_back",object@range,PACKAGE="DuffyNGS") })

setMethod("push_front",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("push_front.bamRange: pushed object must be of class \"bamAlign\"\n")
  .Call("bam_range_push_front",object@range,value@align,PACKAGE="DuffyNGS")
})

setMethod("pop_front",signature="bamRange",definition=function(object)
{.Call("bam_range_pop_front",object@range,PACKAGE="DuffyNGS")})

setMethod("writeCurrentAlign",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("writeCurrentAlign.bamRange: written object must be of class \"bamAlign\"\n")
  .Call("bam_range_write_current_align",object@range,value@align,PACKAGE="DuffyNGS")
})

setMethod("insertPastCurrent",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("insertPastCurrent.bamRange: written object must be of class \"bamAlign\"\n")
  .Call("bam_range_insert_past_curr_align",object@range,value@align,PACKAGE="DuffyNGS")
})

setMethod("insertPreCurrent",signature="bamRange",definition=function(object,value)
{
  if(!is(value,"bamAlign"))
    stop("insertPreCurrent.bamRange: written object must be of class \"bamAlign\"\n")
  .Call("bam_range_insert_pre_curr_align",object@range,value@align,PACKAGE="DuffyNGS")
})

setMethod("moveCurrentAlign",signature="bamRange",definition=function(object,target)
{
  if(!is(target,"bamRange"))
    stop("[moveCurrentAlign.bamRange] target must be bamRange!\n")
  .Call("bam_range_mv_curr_align",object@range,target@range)
  return(invisible())
})


setMethod(f="readID",signature="bamRange",definition=function(object) 
{ .Call("bam_range_get_readID",object@range,PACKAGE="DuffyNGS") })

setMethod(f="readLocus",signature="bamRange",definition=function(object) 
{ .Call("bam_range_get_read_locus",object@range,PACKAGE="DuffyNGS") })

setMethod(f="refID",signature="bamRange",definition=function(object)
{return(.Call("bam_range_get_refid",object@range,PACKAGE="DuffyNGS"))})

setMethod(f="position",signature="bamRange",definition=function(object)
{return(.Call("bam_range_get_position",object@range,PACKAGE="DuffyNGS"))})

setMethod(f="alignLength",signature="bamRange",definition=function(object)
{return(.Call("bam_range_get_align_length",object@range,PACKAGE="DuffyNGS"))})

setMethod(f="alignSeq",signature="bamRange",definition=function(object)
{return(.Call("bam_range_get_align_sequence",object@range,PACKAGE="DuffyNGS"))})

setMethod(f="alignQual",signature="bamRange",definition=function(object)
{return(.Call("bam_range_get_align_qualities",object@range,PACKAGE="DuffyNGS"))})

setMethod(f="readSeq",signature="bamRange",definition=function(object)
{return(.Call("bam_range_get_read_sequence",object@range,PACKAGE="DuffyNGS"))})

setMethod(f="readQual",signature="bamRange",definition=function(object)
{return(.Call("bam_range_get_read_qualities",object@range,PACKAGE="DuffyNGS"))})

setMethod(f="readPhredScores",signature="bamRange",definition=function(object)
{return(.Call("bam_range_get_read_phredScores",object@range,PACKAGE="DuffyNGS"))})

setMethod("getTag",signature="bamRange",definition=function(object, tag)
{ 
  if(!is.character(tag)) stop("[getTag] tag must be character!\n")
  return(.Call("bam_range_get_tag",object@range, tag, PACKAGE="DuffyNGS"))
})

setMethod("getAllTags",signature="bamRange",definition=function(object, sep="::")
{ 
  if(!is.character(sep)) stop("[getAllTags] sep must be character!\n")
  return(.Call("bam_range_get_all_tags",object@range, sep, PACKAGE="DuffyNGS"))
})

setMethod("setTag",signature="bamRange",definition=function(object, tag, value)
{ 
  if(!is.character(tag)) stop("[setTag] tag must be character!\n")
  if ( (n <- size(object)) != length(value)) value <- rep( value, length.out=n);
  return(.Call("bam_range_set_tag",object@range, tag, as.character(value), PACKAGE="DuffyNGS"))
})

setMethod(f="mapQuality",signature="bamRange",definition=function(object)
{ .Call("bam_range_get_map_quality",object@range,PACKAGE="DuffyNGS")})

setMethod("reverseStrand", "bamRange", function(object)
  return(.Call("bam_range_get_strand_reverse",object@range,PACKAGE="DuffyNGS")))

setMethod("paired", "bamRange", function(object)
  return(.Call("bam_range_is_paired",object@range,PACKAGE="DuffyNGS")))

setMethod("firstInPair", "bamRange", function(object)
  return(.Call("bam_range_is_first_in_pair",object@range,PACKAGE="DuffyNGS")))

setMethod("secondInPair", "bamRange", function(object)
  return(.Call("bam_range_is_second_in_pair",object@range,PACKAGE="DuffyNGS")))

setMethod("cigarString", "bamRange", function(object)
  return(.Call("bam_range_get_cigar_string",object@range,PACKAGE="DuffyNGS")))

setMethod("unmapped", "bamRange", function(object)
  return(.Call("bam_range_is_unmapped",object@range,PACKAGE="DuffyNGS")))

setMethod("secondaryAlign", "bamRange", function(object)
  return(.Call("bam_range_is_secondary_align",object@range,PACKAGE="DuffyNGS")))

setMethod(f="mateRefID",signature="bamRange",definition=function(object)
{ .Call("bam_range_get_mate_refid",object@range,PACKAGE="DuffyNGS")})

setMethod(f="matePosition",signature="bamRange",definition=function(object)
{ .Call("bam_range_get_mate_position",object@range,PACKAGE="DuffyNGS")})

setMethod(f="insertSize",signature="bamRange",definition=function(object)
{ .Call("bam_range_get_insert_size",object@range,PACKAGE="DuffyNGS")})

setMethod("mateReverseStrand", "bamRange", function(object)
  return(.Call("bam_range_mate_strand_reverse",object@range,PACKAGE="DuffyNGS")))

setMethod("mateUnmapped", "bamRange", function(object)
  return(.Call("bam_range_mate_is_unmapped",object@range,PACKAGE="DuffyNGS")))

setMethod(f="cigarData",signature="bamRange",definition=function(object, readUnits=FALSE)
{ .Call("bam_range_get_cigar_df",object@range,as.logical(readUnits),PACKAGE="DuffyNGS")})

setMethod(f="mismatchData",signature="bamRange",definition=function(object, readUnits=FALSE)
{ .Call("bam_range_get_mismatch_df",object@range,as.logical(readUnits),PACKAGE="DuffyNGS")})

setMethod("modifyAlign",signature="bamRange", definition=function(object, refID=NULL,
		pos=NULL, mateRefID=NULL, matePos=NULL, readID=NULL, qseq=NULL, qual=NULL)
{ 
  if ( ! is.null( refID)) refID <- as.integer( refID)
  if ( ! is.null( pos)) pos <- as.integer( pos)
  if ( ! is.null( mateRefID)) mateRefID <- as.integer( mateRefID)
  if ( ! is.null( matePos)) matePos <- as.integer( matePos)
  ans <- .Call("bam_range_modify",object@range, refID, pos, mateRefID, matePos, readID, qseq, qual, PACKAGE="DuffyNGS")

  if(is.null(ans))
    	return(ans)
  else
	rangeObject <- bamRange( NULL, NULL, FALSE);
  	rangeObject@range <- ans
	return( rangeObject)
})


###################################################################################################
#                                                                                                 #
# bamAlign                                                                                        #
#                                                                                                 #
###################################################################################################

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# bamAlign encapsulates all contained data in a single dataset in a BAM-file. bamAlign objects    #
# can be read from a bamReader instance and written to a bamWriter instance. All contained data   #
# can be read and written via accessor functions.                                                 #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("bamAlign", representation(align="externalptr"),
         validity=function(object){return(ifelse(is.null(object@align,FALSE,TRUE)))})

setMethod(f="initialize", signature="bamAlign",
          definition=function(.Object,align=NULL){
            .Object@align<-align
            return(.Object)
          }
)


# bamAlign Member Reader functions
setMethod(f="readID",signature="bamAlign",definition=function(object) 
{ .Call("bam_align_get_readID",object@align,PACKAGE="DuffyNGS") })

setMethod(f="readLocus",signature="bamAlign",definition=function(object) 
{ .Call("bam_align_get_read_locus",object@align,PACKAGE="DuffyNGS") })

setMethod(f="refID",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_refid",object@align,PACKAGE="DuffyNGS"))})

setMethod(f="position",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_position",object@align,PACKAGE="DuffyNGS"))})

setMethod("nCigar",signature="bamAlign",definition=function(object)
{ return(.Call("bam_align_get_nCigar",object@align,PACKAGE="DuffyNGS"))})

setMethod(f="cigarData",signature="bamAlign",definition=function(object, readUnits=FALSE)
{ .Call("bam_align_get_cigar_df",object@align,as.logical(readUnits),PACKAGE="DuffyNGS")})

setMethod(f="mismatchData",signature="bamAlign",definition=function(object, readUnits=FALSE)
{ .Call("bam_align_get_mismatch_df",object@align,as.logical(readUnits),PACKAGE="DuffyNGS")})

setMethod(f="mateRefID",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_mate_refid",object@align,PACKAGE="DuffyNGS")})

setMethod(f="matePosition",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_mate_position",object@align,PACKAGE="DuffyNGS")})

setMethod(f="insertSize",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_insert_size",object@align,PACKAGE="DuffyNGS")})

setMethod(f="mapQuality",signature="bamAlign",definition=function(object)
{ .Call("bam_align_get_map_quality",object@align,PACKAGE="DuffyNGS")})

setMethod(f="alignLength",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_align_length",object@align,PACKAGE="DuffyNGS"))})

setMethod(f="alignSeq",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_align_sequence",object@align,PACKAGE="DuffyNGS"))})

setMethod(f="alignQual",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_align_qualities",object@align,PACKAGE="DuffyNGS"))})

setMethod(f="readSeq",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_read_sequence",object@align,PACKAGE="DuffyNGS"))})

setMethod(f="readQual",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_read_qualities",object@align,PACKAGE="DuffyNGS"))})

setMethod(f="readPhredScores",signature="bamAlign",definition=function(object)
{return(.Call("bam_align_get_read_phredScores",object@align,PACKAGE="DuffyNGS"))})

setMethod("getTag",signature="bamAlign",definition=function(object, tag)
{ 
  if(!is.character(tag)) stop("[getTag] tag must be character!\n")
  return(.Call("bam_align_get_tag",object@align, tag, PACKAGE="DuffyNGS"))
})

setMethod("getAllTags",signature="bamAlign",definition=function(object, sep="::")
{ 
  if(!is.character(sep)) stop("[getAllTags] sep must be character!\n")
  return(.Call("bam_align_get_all_tags",object@align, sep, PACKAGE="DuffyNGS"))
})

setMethod("setTag",signature="bamAlign",definition=function(object, tag, value)
{ 
  if(!is.character(tag)) stop("[setTag] tag must be character!\n")
  return(.Call("bam_align_set_tag",object@align, tag, as.character( value), PACKAGE="DuffyNGS"))
})

setMethod("modifyAlign",signature="bamAlign", definition=function(object, refID=NULL,
		pos=NULL, mateRefID=NULL, matePos=NULL, readID=NULL, qseq=NULL, qual=NULL)
{ 
  if ( ! is.null(refID)) refID <- as.integer( refID)
  if ( ! is.null(pos)) pos <- as.integer( pos)
  if ( ! is.null( mateRefID)) mateRefID <- as.integer( mateRefID)
  if ( ! is.null( matePos)) matePos <- as.integer( matePos)
  ans <- .Call("bam_align_modify",object@align, refID, pos, mateRefID, matePos, readID, qseq, qual, PACKAGE="DuffyNGS")
  if(is.null(ans))
    return(ans)
  else
    return(new("bamAlign",ans))
})

setMethod("show","bamAlign",function(object){
  cat("An object of class '",class(object),"'.\n",sep="")
  cat( "readID:   ", readID(object), "\nrefID:    ",refID(object),"\nposition: ",position(object),"\n")
  if( unmapped(object)) {
  	cat( "readSeq: ", readSeq(object),"\nreadQual:", readQual(object), "\n")
  } else {
  	cat( "alignSeq: ", alignSeq(object),"\nalignQual:", alignQual(object), "\n")
  	cat( "mateRefID:    ",mateRefID(object),"\nmatePosition: ",matePosition(object),"\n")
  	cat( "insertSize:", insertSize(obj),"\n")
  	cat("cigarData:\n")
  	print(cigarData(object))
  }
  cat( "tags:    ", unlist( strsplit( getAllTags(object), split="::")))
  cat("\n")
})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Queries against alignment flag (Readers and Accessors)

# pcrORopt_duplicate
setGeneric("pcrORopt_duplicate", function(object) standardGeneric("pcrORopt_duplicate"))
setMethod("pcrORopt_duplicate", "bamAlign", function(object)
  return(.Call("bam_align_is_pcr_or_optical_dup",object@align,PACKAGE="DuffyNGS")))
setGeneric("pcrORopt_duplicate<-", function(object,value) standardGeneric("pcrORopt_duplicate<-"))
setReplaceMethod(f="pcrORopt_duplicate", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, Duplicate setter: value must be boolean")
                   .Call("bam_align_set_is_pcr_or_optical_dup",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# failedQC
setGeneric("failedQC", function(object) standardGeneric("failedQC"))
setMethod("failedQC", "bamAlign", function(object)
  return(.Call("bam_align_fail_qc",object@align,PACKAGE="DuffyNGS")))
setGeneric("failedQC<-", function(object,value) standardGeneric("failedQC<-"))
setReplaceMethod(f="failedQC", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, failedQC setter: value must be boolean")
                   .Call("bam_align_set_fail_qc",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# firstInPair
setMethod("firstInPair", "bamAlign", function(object)
  return(.Call("bam_align_is_first_in_pair",object@align,PACKAGE="DuffyNGS")))
setGeneric("firstInPair<-", function(object,value) standardGeneric("firstInPair<-"))
setReplaceMethod(f="firstInPair", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, FirstInPair setter: value must be boolean")
                   .Call("bam_align_set_is_first_in_pair",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# secondInPair
setMethod("secondInPair", "bamAlign", function(object)
  return(.Call("bam_align_is_second_in_pair",object@align,PACKAGE="DuffyNGS")))
setGeneric("secondInPair<-", function(object,value) standardGeneric("secondInPair<-"))
setReplaceMethod(f="secondInPair", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, secondInPair setter: value must be boolean")
                   .Call("bam_align_set_is_second_in_pair",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)

setMethod("cigarString", "bamAlign", function(object)
  return(.Call("bam_align_get_cigar_string",object@align,PACKAGE="DuffyNGS")))

# unmapped
setMethod("unmapped", "bamAlign", function(object)
  return(.Call("bam_align_is_unmapped",object@align,PACKAGE="DuffyNGS")))
setGeneric("unmapped<-", function(object,value) standardGeneric("unmapped<-"))
setReplaceMethod(f="unmapped", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, unmapped setter: value must be boolean")
                   .Call("bam_align_set_is_unmapped",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# mateUnmapped
setMethod("mateUnmapped", "bamAlign", function(object)
  return(.Call("bam_align_mate_is_unmapped",object@align,PACKAGE="DuffyNGS")))
setGeneric("mateUnmapped<-", function(object,value) standardGeneric("mateUnmapped<-"))
setReplaceMethod(f="mateUnmapped", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, mateUnmapped setter: value must be boolean")
                   .Call("bam_align_set_mate_is_unmapped",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# reverseStrand
setMethod("reverseStrand", "bamAlign", function(object)
  return(.Call("bam_align_get_strand_reverse",object@align,PACKAGE="DuffyNGS")))
setGeneric("reverseStrand<-", function(object,value) standardGeneric("reverseStrand<-"))
setReplaceMethod(f="reverseStrand", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, reverseStrand setter: value must be boolean")
                   .Call("bam_align_set_strand_reverse",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# mateReverseStrand
setMethod("mateReverseStrand", "bamAlign", function(object)
  return(.Call("bam_align_mate_strand_reverse",object@align,PACKAGE="DuffyNGS")))
setGeneric("mateReverseStrand<-", function(object,value) standardGeneric("mateReverseStrand<-"))
setReplaceMethod(f="mateReverseStrand", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, mateReverseStrand setter: value must be boolean")
                   .Call("bam_align_set_mate_strand_reverse",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# paired
setMethod("paired", "bamAlign", function(object)
  return(.Call("bam_align_is_paired",object@align,PACKAGE="DuffyNGS")))
setGeneric("paired<-", function(object,value) standardGeneric("paired<-"))
setReplaceMethod(f="paired", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, paired setter: value must be boolean")
                   .Call("bam_align_set_is_paired",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# properPair
setGeneric("properPair", function(object) standardGeneric("properPair"))
setMethod("properPair", "bamAlign", function(object)
  return(.Call("bam_align_mapped_in_proper_pair",object@align,PACKAGE="DuffyNGS")))
setGeneric("properPair<-", function(object,value) standardGeneric("properPair<-"))
setReplaceMethod(f="properPair", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, properPair setter: value must be boolean")
                   .Call("bam_align_set_mapped_in_proper_pair",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# secondaryAlign
setMethod("secondaryAlign", "bamAlign", function(object)
  return(.Call("bam_align_is_secondary_align",object@align,PACKAGE="DuffyNGS")))
setGeneric("secondaryAlign<-", function(object,value) standardGeneric("secondaryAlign<-"))
setReplaceMethod(f="secondaryAlign", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.logical(value))
                     stop("class bamReader, SecondaryAlign setter: value must be boolean")
                   .Call("bam_align_set_is_secondary_align",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
# flag
setGeneric("flag", function(object) standardGeneric("flag"))
setMethod("flag", "bamAlign", function(object)
  return(.Call("bam_align_get_flag",object@align,PACKAGE="DuffyNGS")))
setGeneric("flag<-", function(object,value) standardGeneric("flag<-"))
setReplaceMethod(f="flag", signature="bamAlign",
                 definition=function(object,value){
                   if(!is.integer(value))
                     stop("class bamReader, flag setter: value must be boolean")
                   .Call("bam_align_set_flag",object@align,value,PACKAGE="DuffyNGS")
                   return(object)
                 }
)
#  End: Queries against alignment flag (Readers and Accessors)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

#  End: bamAlign
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  coercing

as.data.frame.bamRange<-function(x,row.names=NULL,optional=FALSE,...)
  {return(.Call("bam_range_get_align_df",x@range,PACKAGE="DuffyNGS"))}

as.data.frame.refSeqDict<-function(x,row.names=NULL,optional=FALSE,...)
{
  if(is.null(row.names))
    row.names<-1:(length(x@SN))
  else if(length(row.names)!=length(x@SN))
    stop("[as.data.frame.refSeqDict] length(row.names)!=length(x@SN)!")
  return(data.frame(SN=x@SN,LN=x@LN,AS=x@AS,M5=x@M5,SP=x@SP,UR=x@UR,row.names=row.names))  
}


setAs("bamRange","data.frame",function(from)
  {return(.Call("bam_range_get_align_df",from@range,PACKAGE="DuffyNGS"))})
setAs("refSeqDict","data.frame",function(from)
  {return(data.frame(SN=from@SN,LN=from@LN,AS=from@AS,M5=from@M5,SP=from@SP,UR=from@UR,row.names=1:length(from@SN)))})



# +++++++++++++++++++++++++++++
# helper functions

refID2seqID <- function( refID, reader=NULL, refData=NULL) {

	# 2 ways to get to the ref data table
	if ( is.null(refData) && ( ! is.null(reader))) {
		refData <- getRefData( reader)
	}
	if ( is.null(refData)) stop( "'refID2seqID' needs non-null reader or refData object")

	out <- rep( NA, times=length(refID))
	hits <- match( refID, refData$ID, nomatch=0)
	out[ hits > 0] <- refData$SN[hits]

	return( out)
}


seqID2refID <- function( seqID, reader=NULL, refData=NULL) {

	# 2 ways to get to the ref data table
	if ( is.null(refData) && ( ! is.null(reader))) {
		refData <- getRefData( reader)
	}
	if ( is.null(refData)) stop( "'seqID2refID' needs non-null reader or refData object")

	out <- rep( -1, times=length(seqID))
	hits <- match( seqID, refData$SN, nomatch=0)
	out[ hits > 0] <- refData$ID[hits]

	return( out)
}


`softClipAlignLength` <- function( alignlens, cigars) {

	# if there was any soft clipping, there are 'S' in the cigar strings
	whoClip <- grep( "S", cigars, fixed=T)
	if ( ! length(whoClip)) return( alignlens)

	out <- alignlens

	clipLeft <- grep( "^[0-9]+S", cigars[whoClip])
	if ( length(clipLeft)) {
		whoLeft <- whoClip[ clipLeft]
		digits <- sub( "(^[0-9]+)(S.+)", "\\1", cigars[whoLeft])
		nClipLeft <- as.numeric( digits)
		out[ whoLeft] <- out[ whoLeft] - nClipLeft
	}
	clipRight <- grep( "[0-9]+S$", cigars[whoClip])
	if ( length(clipRight)) {
		whoRight <- whoClip[ clipRight]
		digits <- sub( "(^.+[MIDNPX=])([0-9]+)(S$)", "\\2", cigars[whoRight])
		nClipRight <- as.numeric( digits)
		out[ whoRight] <- out[ whoRight] - nClipRight
	}
	return(out)
}
