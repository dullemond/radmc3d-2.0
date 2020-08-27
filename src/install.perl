#!/usr/bin/perl
#=======================================================
#           INSTALL ROUTINE FOR "RADMC-3D"
#=======================================================
$pwd  = `pwd` ; 
chop($pwd) ;
print "Installing RADMC-3D in the directory:\n  '$pwd'\n";
$home = $ENV{"HOME"};
$bin  = $home . "/bin" ;
if(!(-e $bin)) {
    print "You must have a bin/ directory in your home directory\n" ;
    print "-----Shall I make $bin for you?\n" ;
    $input = <STDIN> ;
    print $input ;
    if($input=~/^[yY]/) {
	print "Creating $bin for you...\n" ;
	system("mkdir $bin") ;
    }
}
$path = $ENV{"PATH"};
if(!($path =~ /$bin/)) {
    print "The $bin directory exists, but it not in the PATH environment variable\n" ;
    print "You must put the \n" ;
    print "$bin directory \n" ;
    print "in the path yourself (in the .tcshrc file if you use the tcsh shell...)\n" ;
    print "If you do it now, don't forget to type 'rehash'.\n" ;
}
print "  Creating a link 'radmc3d' in '$bin/'\n" ;
$radmc3d    = $pwd . "/radmc3d" ;
$radmc3dlnk = $bin . "/radmc3d" ;
if(!(-e $radmc3dlnk)) {
    print "------ Warning: file $radmc3dlnk did not exist previously. You might want to type 'rehash'\n" ;
}
open(FILE,">$radmc3dlnk") || die "Could not open file\n" ;
print FILE "#!/usr/bin/perl\n" ;
print FILE "system(\"$radmc3d \@ARGV\");" ;
close (FILE) ;
`chmod u+rwx $radmc3dlnk` ;
#
# The following is only necessary if you want to use the auxiliary Python 
# tools in the python/radmc3d_tools/ directory. You can remove or comment-out the 
# following lines if you do not want to use these tools.
#
# NOTE: For the radmc3dPy package (in the python/radmc3dPy/ directory), 
#       you should use the installation method of that package. Read the
#       python/radmc3dPy/README file for directions.
#
$python    = $home . "/bin/python" ;
$radmc3d   = $python . "/radmc3d_tools" ;
if(!(-e $python)) {
    print "It is convenient for you to have a bin/python/ directory in your home directory\n" ;
    print "-----Shall I make $python for you?\n" ;
    $input = <STDIN> ;
    print $input ;
    if($input=~/^[yY]/) {
	print "Creating $python for you...\n" ;
	system("mkdir $python") ;
    }
}
if(-e $python) {
    if(!(-e $radmc3d)) {
        print "It is convenient for you to have a bin/python/radmc3d_tools directory in your home directory\n" ;
        print "-----Shall I make $radmc3d for you?\n" ;
        $input = <STDIN> ;
        print $input ;
        if($input=~/^[yY]/) {
            print "Creating $radmc3d for you...\n" ;
            system("mkdir $radmc3d") ;
        }
    }
    if(-e $radmc3d) {
        print "Found the ~/bin/python/radmc3d_tools directory. Copying current python radmc3d tools there.\n" ;
        system("rm -f $radmc3d/*.pyc") ;
        system("cp -r ../python/radmc3d_tools/* $radmc3d/") ;
    }
}
