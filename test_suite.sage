attach("mf_secret_sharing.sage")

def test_recon_random(level, weight, prime, num_senders, r):
    ring = GF(prime)
    M = ModularForms(level, weight, base_ring = ring)
    pub = create_public_data(M, num_senders, r)
    votes = [ring(randint(0,1)) for ix in range(0,num_senders)]
    sum_votes = sum(votes)

    #senders send (partial) votes to auditors
    auditors = [ring.zero() for ix in range(0,M.sturm_bound()+1)]
    for sender in range(0,num_senders):
        for auditor in range(0,len(auditors)):
            #auditors collect (partial) votes
            auditors[auditor] += send_to_auditor(votes[sender], sender, auditor, pub)

    #auditors send (partial) votes to senders
    sender_recv=[[] for ix in range(0,num_senders)]
    for auditor in range(0,len(auditors)):
        for sender in range(0,num_senders):
            #need id because order is important, and comms might be async
            sender_recv[sender].append((auditor,auditors[auditor]))

    #senders reconstruct vote sum
    sum_votes_recon = [ring.zero() for ix in range(0,num_senders)]
    for sender in range(0,num_senders):
        sum_votes_recon[sender] = reconstruct(sender, sender_recv[sender], pub)
        assert sum_votes_recon[sender] == sum_votes, "Reconstructed vote sum is incorrect."

def runtests_random_r():
    for (lvl,wgt,pr) in [(23,2,7),(23,2,11),(29,2,7),(29,2,11),(127,2,7),(127,2,11),(127,2,13)]:
        M = ModularForms(lvl,wgt,base_ring=GF(pr))
        for num_senders in range(1,min(pr-1,M.dimension())+1):
            r = randint(100,200)
            try:
                test_recon_random(lvl,wgt,pr,num_senders,r)
            except ValueError as err:
                print("ValueError; {0}".format(err))
            except RuntimeError as err:
                #this is BAD!
                print("RuntimeError; {0}".format(err))
