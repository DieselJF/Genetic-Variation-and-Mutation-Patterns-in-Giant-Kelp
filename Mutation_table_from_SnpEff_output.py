# This short script will parse the filtered, decomposed, and SnpEff annotated VCF file to a more human readable table.
# The final mutation table columns are
# 1- Gene: Name of affected gene
# 2- Scaffold: Scaffold where mutation is found
# 3- Position: Starting nucleotide affected by mutation
# 4- Ref: Reference allele
# 5- Alt: Alternate allele
# 6- Effect: Effect of the variant
# 7- Individuals: Individuals in the seedbank containing the mutation
# 8- SnpEff annotation: Full SnpEff annotation containing effect, impact, functional class, codon change, amino acid change, etc. (more here: http://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files)

# Opening vcf file
vcf = open("snpEff_annotated_decomposed_haploid_559_indv_on_CI_03_genome_final_indel_and_snp_gatk4_hard_filter_minq30_minmeanDP2.vcf", "r")

#Index based on vcf input file
indexes = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1	10	101	102	103	104	105	106	107	108	109	11	110	111	112	113	114	115	116	118	119	12	120	121	122	123	124	125	126	127	128	129	13	130	131	132	133	134	135	136	137	138	139	14	140	141	142	143	144	145	146	147	148	149	15	150	151	152	153	154	155	156	157	158	159	16	160	161	162	163	164	165	166	167	168	169	17	170	171	172	173	174	175	176	177	178	179	18	180	181	182	183	184	185	186	187	188	189	19	190	191	192	193	194	195	197	198	199	2	20	200	201	202	203	204	205	206	207	208	209	21	210	211	212	213	214	215	216	217	218	219	22	220	221	222	223	224	225	226	227	228	229	23	230	231	232	233	234	235	236	237	238	239	24	240	241	242	243	244	245	246	247	248	249	25	250	251	252	253	254	255	256	257	258	259	26	260	261	263	264	265	266	267	268	269	27	270	271	272	273	274	275	276	277	278	279	28	280	281	282	283	284	285	286	287	288	289	29	290	291	292	293	294	295	296	297	298	299	3	30	300	301	302	303	304	305	306	307	308	309	31	310	311	312	313	314	315	316	317	318	319	32	320	321	322	323	324	325	326	327	328	329	33	330	331	332	333	334	335	336	337	338	339	34	340	341	342	343	345	346	347	348	349	35	350	351	352	353	354	355	356	357	358	359	36	360	361	362	363	364	365	366	367	368	369	37	370	371	372	373	374	375	376	377	378	379	38	380	381	382	383	384	385	386	387	388	389	39	390	391	392	393	394	395	396	397	398	399	4	40	400	401	402	403	404	405	406	408	409	41	410	411	412	413	414	415	416	417	419	42	420	421	422	423	424	425	426	427	428	429	43	430	431	432	433	434	435	436	437	438	439	44	440	441	442	443	444	445	446	447	448	449	45	450	451	452	453	454	455	456	457	458	459	46	460	461	462	463	464	465	466	467	468	469	47	470	471	472	473	474	475	476	477	478	479	48	480	481	482	483	484	485	486	487	488	489	49	490	491	492	493	494	495	496	497	498	499	5	50	500	501	502	503	504	505	506	507	508	509	51	510	511	512	513	514	515	516	517	518	519	52	520	524	526	527	528	529	53	531	532	533	534	535	536	537	539	54	540	541	542	543	544	545	546	547	548	549	55	550	551	552	553	554	556	557	558	559	56	561	562	563	564	565	566	567	57	570	571	572	573	574	575	576	577	578	579	58	580	581	586	6	60	62	63	64	65	66	67	68	69	7	70	71	72	74	75	76	77	78	79	8	80	82	83	84	86	87	88	9	90	91	92	93	94	95	96	97	98	99"
indexes = indexes.split()

fileOut = open("Mutation_results_snpEff_annotated_decomposed_haploid_559_indv_on_CI_03_genome_final_indel_and_snp_gatk4_hard_filter_minq30_minmeanDP2.txt", "w")

for line in vcf:
    if not line.startswith("#"):
        line = line.strip().split()
        scaffold = line[0]
        pos = line[1]
        Ref_al = line[3]
        Alt_al = line[4]
        if Alt_al == "*":
            continue
        else:
            gene_eff_list = line[7].split(";") #
            for item in gene_eff_list:
                if item.startswith("ANN"):
                    item = item[4:]
                    gene_eff_list = item.split(',')

            for item in gene_eff_list:
                if item.startswith("("): # Correct LOF, (comma) that shows up in just a few cases
                    continue
                else:
                    info = item
                    effect = item.split("|")[1]
                    gene = item.split("|")[4] + "|" + item.split("|")[5]
                    individuals_index = [i for i in range(len(line)) if line[i].startswith("1:")]
                    individuals_ref = [indexes[i] for i in individuals_index]
                print(gene, scaffold, pos, Ref_al, Alt_al, effect, ",".join(individuals_ref), info, sep="\t") #, file=fileOut)
vcf.close()


