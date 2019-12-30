# ubiome数据分析流程学习笔记1

从科研的角度讲，肠道微生物的研究依然大热，cns大作文章层出不穷，带来新的idea和见解，另一方面，微生物产业却道路曲折，根据肠道产业公众号的报道：

> **uBiome 仅以 1%残值出售知识产权**
>
> 2019 年 12 月 17 日，已经在 9 月 4 日向美国破产法院申请破产保护的 uBiome 公司宣布，将其名下 246 项微生物相关专利和 30 万人微生物组数据，以 770 万美元的价格转让给 Psomagen 公司，后者是韩国生物技术公司 Macrogen 在美国开展基因检测业务的子公司。

一般Next-Seq等illumina机器的最具性价比的测序读长是2＊150（150PE），而16S-V4区的长度是806R－515F＝291bp。那么，使用150PE测序方法对16S-V4区测序的话，获得的测序数据重叠区的长度是150＋150－291（包括引物长度）＝9bp，这个重叠理论上是可以拼接的，可是实际上，由于illumina测序反应（每次SBS会延伸一个碱基，大约耗时70分钟，据说之前拍照就需要40分钟。2013，https://blog.csdn.net/seallama/article/details/17613601） 推测大概在几十分钟，所以，150＊2个反应下来，试剂的衰减相当严重，测序最后获得的数据质量也急剧下降了，想要通过这最多9个的重叠碱基来实现拼接几乎是不可能完成的任务。

可是，这基本上是许多测序公司的测序方法，如果他们不采用特殊的方法是不能完成拼接，或者实现双向数据的有效利用的。其实这个问题，有两个方向可以解决，都已经发表了论文，一个是使用16S－V4区的通用引物来取代illumina的p5，p7测序引物，这样有效测序长度就变为300bp左右（不算上正反20bp左右的引物），重叠区域就变为50bp，许多论文和公司就是这么做的，这样做也比较简单和可靠，对于分析难度也比较小。如果一个lane全部测16S-V4，这是没问题的。

然而，16S-V4的测序数据量一般比较小，即使对于通量比较小的NextSeq－500机器，比如每次产出20GB算，一个样本获得测序rawdata 50Mb（量应该足够多了）算的话，至少需要500个样本混样测序（可能为了增加序列多样性，可能还要增加一定比例的phi噬菌体序列），而这对一般的检测机构，这么大的量是不可思议的。如果和别的不是测16S-V4文库一起测序的话就是不可行的，特别是16S－V4文库占比比较小的情况下。

这里，就引出了另一个方法，方法不好实现，就只有技术来扛了。另一个解决方案便是通过高水平的分析来解决不能拼接的问题了。我找到了几个方法，其中之一是读到的ubiome公司的一个方法。鉴于ubiome公司已经破产，据说是由于医保骗保，还有小道消息说是FBI盯上了遗传资源，但是从上面的报道来看，数据都可以卖给韩国人，小道消息应该是不正确。这个方法仅供参考，而且，由于他们公司注册了相关专利，也只能是学习下了，不可能商业应用。这篇文章发表在plos one，并不是多么高级的期刊，而且，这个公司的检测结果也面临质疑，所以，还是仅供参考。

## **数据分析流程**

（其申请了专利，流程较复杂，特别是数据库要结合实验处理以减少假阳性和假阴性。）

### **1)数据库准备**

a.首先从SILVA-16S数据库中找出能用V4通用引物扩出的序列，允许两个错配。检查得到的序列里是否有简并碱基，>20个以上的序列质量较差，由于采用双端测序，模拟测序产生的序列，去掉引物，就相当于正向有125bp和反向124bp接在一起。

b. 通过有选择地删除每个分类群的非特定扩增子集，可以为每个分类群创建几个经过筛选的数据库，并使用下面概述的步骤确定最佳数据库。采用100%的序列相似度和长度进行分析，排除不特异的扩增序列。然后，将注释到感兴趣分类群(Ti)和不同分类群(Dt)的序列按照dt/ti排列，并根据它们的商值分组。例如，在一个极端上，当商数为0时，我们创建了非常具体的数据库，其中不允许出现假阳性。在另一个极端，当商数很大时，例如，100，我们创建了非常敏感的数据库，在那里我们最大限度地识别真阳性，代价是产生假阳性。我们广泛地探索了这两个极端之间的不同可能的数据库，使用每个生成的数据库来预测每个感兴趣的物种和属的所有16S序列的分类。评佑它们的性能，换句话说，选择性地从数据库中删除带有非特异性扩增的序列，同时最大限度地提高识别每个分类组中大多数序列的敏感性、特异性、精密度和阴性预测值。选择了灵敏度、特异度、精密度和阴性预测值均在90%以上、精密度和特异度之间的距离为最小可能值的分类单元作为每个分类单元的最佳数据库，目的是尽可能使精密度优于特异度。最终选择了28个种属进行分析。

### **2）物种注释**

a.Q30过滤

b.去除引物和接头等

c. Swarm 2.15聚类，参数使用1个核苷酸的距离，fastidious和usearch-abundance。每个簇中最丰富的序列被认为是真正的生物序列，并被指定为簇中所有读取的计数。聚类中的其余读取被认为包含错误的测序的结果。

d. 使用vSearch从所有聚类中删除嵌合体。

e. 使用100%长度和100%同一性与来自手工处理过的silva数据库123版目标16Srna基因序列和分类注释的数据库进行比对。

f. 每个分类群的相对丰度是通过将与该分类群相关联的计数除以过滤后的总reads来确定的。

## **3.实验验证**

从xTAG GPP (Luminex‘s xTAG Gastrointestinal Pathogen Panel)获得的样品, 对生物信息学流程进行了优化，以确保阳性结果真正意味着目标存在于样本中，并且仅在样本中没有目标存在时才获得阴性结果。

## **4.** **健康队列的参考范围**

为了确定28个种属的健康参考范围，我们建立了897健康个样本的队列。对于897个样本中的每一个，测定了微生物种群中每个目标的相对丰度。这些数据用来定义一个健康参考范围，例如艰难梭状芽孢杆菌在健康队列的~2%中发现，因此我们为其定义了一个相对丰度为0%到0.18%的健康范围。

## **5.不足之处**

该方法不能识别肠球菌内不同的血清型，也不能检测能够区分致病性艰难梭菌或大肠杆菌与非致病性菌株的毒素基因，也不能在某些属水平的靶标中分辨物种。

## **6.** **报告里的菌的情况**

以上内容为文章中的，示例报告里的内容比文章中的有更新，菌有64属（种），35个种，29属。

可以用snp芯片数据实现HLA分型，任意的SNP芯片，只要位点数足够即可（几十万）。由于本人几乎没有前端和后
[qiime2+picrust1学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380311&idx=2&sn=c96b26c9085f70a7b683c29e5d90c6e3&chksm=817095fab6071cec75ac28103042fa570285c852fc99220ec9d8549f5c4e62777f935ab3297f#rd)
[画个草原之旅路线图](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380303&idx=1&sn=6026f47e2ccda39bd8e1dd935bea8326&chksm=817095e2b6071cf41bdbb1b680122d1d1a7fa7ca8f204d3c651fffbb7b5abb00b0a96f060428#rd)
[欣赏下内蒙风光](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380303&idx=2&sn=2c3382257652cd85ac360d72c446427b&chksm=817095e2b6071cf4e23a74809fcd7f2489d5324d714375460f471a531796bd87a3f2b9a99799#rd)
[纳米孔Nanopore-16S数据分析学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380285&idx=1&sn=d5525dae50e23ef3bcb38a2ec8f3ea8d&chksm=81709590b6071c8697c0f3606db006f3acfd677100c38fc422ff33fd2e4cc814824408f924aa#rd)
[Nanopore牛津纳米孔测16S学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380281&idx=1&sn=dbc9ba5820c0e32bc820902fa637db87&chksm=81709594b6071c8208cf24ea5ee59e9dfb8990dee92e0aba48b611fa06fcf3430da41d4f8e19#rd)
[2019年主流测序仪参数对比](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380267&idx=1&sn=589f3af8d5d8574acf1c4aeb926589fa&chksm=81709586b6071c90e4460e4ce11dc052671b6d5710a932c816a339d1756f01aa30230e15b05a#rd)
[qiime2-2019.4更新学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380260&idx=1&sn=3eb486ea668638b970fed8b671cd8a11&chksm=81709589b6071c9f795a2c624770d91644c081979bd8a4e9e9e64f4e3bcc23e9b32d10b746a6#rd)
[祖源探索"三步曲"](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380254&idx=1&sn=f5ca938b5a37b1398f12dfffe6a60430&chksm=817095b3b6071ca5cabb3d636d63558f6340ed12e8fc4d160781e7971d70e5f80606ff1e7fbb#rd)
[体验impute.me基因检测分析结果](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380249&idx=1&sn=b271d4e485e2b2e34106589931dfbf68&chksm=817095b4b6071ca2a05d112f766a3b89eca3e2918d1a770fa90e5d523b67c639ba3828bfa151#rd)
[qiime2-2019.1更新学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380245&idx=1&sn=5243c586c7e4cdeeb74feb5deee1f2c1&chksm=817095b8b6071cae82cac2fb4d739addefc9f54ff01baf92cdb231ad28a999c152ffa157e9b8#rd)
[肠型分析学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380241&idx=1&sn=fc3bf696c9af066ec66fe440de9638dd&chksm=817095bcb6071caadac43575bef9dcd69290e134c2085de6b426250e97a0fe37f8463b4bc22e#rd)
[qiime2-2018.11发布学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380236&idx=1&sn=6af7c92bb51667c44eeb8a03141396d9&chksm=817095a1b6071cb721e3277848eb70fc323dd75f8e9805f1d2f0ea749f7f680d1e5f53152d0e#rd)
[《让身体和微生物成为朋友》读书笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380231&idx=1&sn=f2c7eeba47739b879bec619bbfd6f87b&chksm=817095aab6071cbc8185dad8b5172bad5504386f7329c9e7da4be8c11c9b50bb05f320c638d5#rd)
[QIIME2-CLI更新学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380226&idx=1&sn=12dcba530c40d5ce021e46b93dca21a5&chksm=817095afb6071cb9dc5cb42fbb4c2c0a9f23c82cc2bec86517d2bdaa66fd0e81866272633aae#rd)
[初识OXFORD NANOPORE](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380223&idx=1&sn=1334fb4e924140a00263db365668bf59&chksm=81709552b6071c4404ab745c98e9c5b7c7be297e53f269db8d09f74760614ef990b66c61a2bb#rd)
[GUTHEALTHY镜像使用笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380219&idx=1&sn=bb327f3d3ad5202c01e5e36fa0b582df&chksm=81709556b6071c40887decaaeec8c1a7892e0ca1441dc52d7a1e03f28c486a6c2eff5ca0f188#rd)
[SILVA、GREENGENES、RDP三大数据库的序列探索统计](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380216&idx=1&sn=4f15fc8741fbcf7114466df97cf64c3a&chksm=81709555b6071c43f129acf515a1980579fc43619befa32939a6964228472356b9c2d1242062#rd)
[做了个简陋的网页娱乐版HLA分型](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380212&idx=1&sn=4e7f58ac60aeeff96d5adc93226f0ae9&chksm=81709559b6071c4f40eb6ff9bf71505e477d8bf7b758707e233373152de6c1049752ed30c4ed#rd)
[SNP2HLA学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380208&idx=1&sn=2d62a300e0592e7c39e219c7cb20d6ce&chksm=8170955db6071c4b7069df6ece33a1a1fe1591e8279884156afb8a7b5a977ff2f828044d3ec5#rd)
[探索一个消费级基因检测结果](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380203&idx=1&sn=b994746ea29b3c11083e96595d660111&chksm=81709546b6071c501840b12c829ad4d2381b1d9339fbc0dda0337d23c098b86f205fc2763061#rd)
[几种有趣的益生菌](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380198&idx=1&sn=409c3aa9b9a4e84c6366b959cc43aca5&chksm=8170954bb6071c5d05bd9827e89b2129f49973eb0c7243121a3550a89c6d75562280c84280b4#rd)
[肠道菌群学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380195&idx=1&sn=f1a52bb5ab5019bdd0f4fd396433db5c&chksm=8170954eb6071c581e03337c65806a0c2b442728e4073e566fd4acdc936dc5778741359a22f3#rd)
[几个肠道微生物检测方法的对比](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380192&idx=1&sn=b11782e28f2fda6c22982210468e839f&chksm=8170954db6071c5b2fdbbedae6ea1a58954e2c6ec6ef7d0c9c8c9231ac3efe7838b06e6735fc#rd)
[Q2studio学习笔记三](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380185&idx=1&sn=8478781403c0a95cd346dca1c235211d&chksm=81709574b6071c624630dbb08a1182d58c2eea6557735f842a9796ae7fbc5e6b181439b9c4c4#rd)
[QIIME2图形界面版（Q2STUDIO）学习笔记二](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380182&idx=1&sn=bb202242263b4b9d327820de94103af0&chksm=8170957bb6071c6d9722082762346683d0dfeb8b02bdff416483870d2057100c20286c456969#rd)
[qiime2图形界面安装学习笔记](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380179&idx=1&sn=04904bf08dae6092c304ecabe7108bb3&chksm=8170957eb6071c68737a3db07368e7b06f66df5938aef3efdbc97fb09aeaf0e1bcf18f49013f#rd)
[QIIME2学习笔记（二）](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380173&idx=1&sn=fbea710a608f9d5ceb9092663069fd32&chksm=81709560b6071c76b8933a7316d4d23ce1386c6a4b8e6907e5d0fc3fa2efd6fbf951d1dc89c7#rd)
[qiime2学习笔记（一）](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380169&idx=1&sn=d06e0e491cc1e8b45aafcdafc6c5b114&chksm=81709564b6071c72bb65d7adc5cee6f736cc8eb7d2661fab0d1b3793ad688e89950d7b914f81#rd)
[使用R语言读取PUBMED存入MYSQL数据库](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380162&idx=1&sn=de5ee9abc42673a0dbe792c09d19edd1&chksm=8170956fb6071c7989c1992b3618a35b3b4eef702e07e62e2e93cf4c1eb5884e6274ca0afc2d#rd)
[翻译--肺微生物组与肺癌之间的相互作用](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380158&idx=1&sn=0f1c326b6e2813300583168e116a3e43&chksm=81709513b6071c0590ca199721c4bea479c4478e12a465cce068a5dd90690e2f41b89551a492#rd)
[微生物群对肺内稳态和疾病的作用](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380153&idx=1&sn=9b1137d06a21487ff244b26b6ba4c569&chksm=81709514b6071c025c895cd1c6177883d7b6933932c65e928717a8988795174feea96039fd75#rd)
[环境微生物学课题组揭示除草剂麦草畏微生物代谢机制及其应用前景](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380149&idx=1&sn=7e6015f2b679548c56842b17efecccb8&chksm=81709518b6071c0e56f382ae52e99acfbbf75e8e16724c9b9508b17f0bf6cac2144ef40edca9#rd)
[kindle dx当显示器，凑合着用吧](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380140&idx=1&sn=2c06143d75d73f0c4665e745fae7c4e7&chksm=81709501b6071c1707125885498a1855becda03328581b0087d39b8f78dde66a8c1c83ce72c6#rd)
[3个肠道微生物检测报告的对比--你有兴趣做肠道微生物检测吗？](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380137&idx=1&sn=2c59b8cdfe29bcdc5348953f6619b7e0&chksm=81709504b6071c128850b132f60d2da5a257f22533b3684bed8704c174731abf98796c39125d#rd)
[Roche480软件安装，win10亲测可用](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380133&idx=1&sn=76998998d34f4ee59238476126d01e59&chksm=81709508b6071c1efeec99d3a7bcbe8042c6c81b05559c29de40e60e0ecb4679a8dcacbc4c04#rd)
[HLA NGS数据分析](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380130&idx=1&sn=70e8cd1d21158cc127fd6271c6600aa1&chksm=8170950fb6071c19db29a163310e7f7d7f94b2843575713251407e4379fcb1190471ce7e016e#rd)
[高分辨率系统发育微生物群落剖析](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380127&idx=1&sn=a61c27195016f3bc8a4c455f4c94ea77&chksm=81709532b6071c248e13656e1f54d460a73fe20f772a9cea226049f550960f536e9df4b4c300#rd)
[分析粪便微生物移植后患者高通量单分子实时测序数据的工作流程](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380127&idx=2&sn=7b351b39301399cdbc2b3edfdd4d5c78&chksm=81709532b6071c247f8912b359c14231fa9e1681dec12a79535f01ba6edf0713088b1848f03a#rd)
[使用CCS序列数据改进宏基因组拼接效率和物种分类注释](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380127&idx=3&sn=aa3dc83ef2213a6a1fd8f0f7a33075d6&chksm=81709532b6071c2475eac0ce0437a8d7c37630b89727c33846e14ac4b2b560ff732ec6200dff#rd)
[科技记者们的门派](http://mp.weixin.qq.com/s?__biz=MzIwMDQ3Njk5NA==&mid=2457380127&idx=4&sn=fbd25f44f2eb1fc030e40df850a86cbf&chksm=81709532b6071c244fd27b07a9c1dfe8e014135074bd85235daa4adc09170323a90694382a09#rd)