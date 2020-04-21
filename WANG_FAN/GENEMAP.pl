#! /usr/bin/perl -w  
use strict;  
use DATA::Dumper;

open (ft_gpl,"GSE19804_family.soft")||die "open error1!\n";
open (ft_rma,"GSE19804_rma.txt")||die "open error2!\n";
open (f_f,">GSE19804.txt")||die "open error3!\n"; 

my @tmp;
my @gene;
my $i;
my %h_gpl;
my %h_symbol;
my %hash_count;
my $sym;
my $symbol;
my $len;

chomp(my $header1=<ft_gpl>);
chomp(my $header2=<ft_rma>);
my @header=split(/\t/,$header2);
my $sample=@header;
print $sample;


while (<ft_gpl>){#逐行读取文件内容
	chomp;
	if(!/^#(.*)/ && !/^!(.*)/){
		@tmp=split(/\t/,$_);	
		#print f_f1 "tmp\t",$tmp[0],"\t",$tmp[10],"\n";
		if(exists ($tmp[10]) ){  ####如果$tmp[10]存在，即基因symbol存在时
			$h_gpl{$tmp[0]}=$tmp[10];
		}
	}
}
#print Dumper(\%h_gpl);

while(<ft_rma>){ #若有相同探针，则取其表达值均值
	chomp;
	@tmp=split(/\t/,$_);
	$tmp[0]=~s/\"//g;
	#print @tmp[0],"\n";
	if( exists ($h_gpl{$tmp[0]})){
		$symbol=$h_gpl{$tmp[0]};
		if($symbol){
			if( exists $h_symbol{$symbol}){# %h_symbol以$symbol为key,表达值为value[如果symbol有相同探针时]
				for($i=1;$i<=$#tmp;$i++){
					$h_symbol{$symbol}->[$i]+=$tmp[$i];
				}
				$hash_count{$symbol}+=1;	# %hash_count以$symbol为key,探针数为value(symbol有相同探针，对探针数计数)	
			}
			else{ #[如果symbol没有相同探针时]
				for($i=1;$i<=$#tmp;$i++){
					$h_symbol{$symbol}->[$i]=$tmp[$i];
				}
				$hash_count{$symbol}=1;	
			}
		}
	}	
}

for($i=1;$i<$sample;$i++){
	print f_f "\t$header[$i]";  #加入表头
}
print f_f "\n";

foreach $sym(keys %h_symbol){
	$len=keys %h_symbol;
	if(!($sym=~/\/\/\//)){   #去///后面内容 
		print f_f $sym;
		for($i=1;$i<$sample;$i++){
			$h_symbol{$sym}->[$i]/=$hash_count{$sym};	
			print f_f "\t$h_symbol{$sym}->[$i]";
		}
		print f_f "\n";
	}
	

}
	
close (ft_gpl);  #关闭数据源文件
close (ft_rma);  
close (f_f);  