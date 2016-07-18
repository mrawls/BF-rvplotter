c
c=======================================================================
c
c  string appendaging
c  the difference between this routine and routine fappend is the "."
c
      subroutine          sappend(instring,appstring,outstring)   
c
c
c
c-----------------------------------------------------------------------
c
c  append appstring on instring and output the result in outstring
c
c=======================================================================
      integer             i,k,lend      
      character*(*)       instring,appstring,outstring
c
      if (appstring.eq.' ') then
       write(6,*) 
     @      ' -WARNING(sappend)-: appstring is blank: no action'
       return
      end if
      if (instring.eq.' ') then
       write(6,*) 
     @      ' -WARNING(sappend)-: instring is blank: no action'
       return
      end if
c
       lend = len(instring)
       k    = 0
       do 09 i=1,lend
        k = i
        if (instring(i:i).eq.' ') goto 10
 09   continue
c
 10   outstring = instring(1:k-1)//appstring      
c
c  normal exit
       return
c
       end
c
c=======================================================================
c
c  file name suffix replacement or appendaging
c  the difference between this routine and routine sappend is the "."
c
      subroutine          fappend(infile,delim,outfile)   
c
c
c
c-----------------------------------------------------------------------
c
c  replace the delimeter on infile with delim and output the
c  result in outfile
c  if the infile has no delimeter, then append delim prefixed 
c  with a "." to infile and output the result in outfile
c
c=======================================================================
      integer             i,k,lend      
      character*(*)       infile,delim,outfile
c
      if (infile.eq.' ') then
       write(6,*) 
     @      ' -WARNING(fappend)-: instring is blank: no action'
       return
      end if
c
       lend = len(infile)
       k    = 0
       do 09 i=1,lend
        k = i
        if ((infile(i:i).eq.'.').or.(infile(i:i).eq.' ')) goto 10
 09   continue
c
 10   if (delim.eq.' ') then
       outfile = infile(1:k-1)
      else
       if (infile(k:k).eq.'.') outfile = infile(1:k)//delim
       if (infile(k:k).eq.' ') outfile = infile(1:k-1)//"."//delim
      end if
c
c
c  normal exit
       return
c
       end
c
c=======================================================================
